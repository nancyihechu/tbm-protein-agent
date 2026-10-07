"""Metric/correspondence tests using analytic geometries and local PDB fixtures."""
from __future__ import annotations

from io import StringIO
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from Bio.PDB import MMCIFIO, PDBParser

from agent2.metrics import (MappingError, _write_ca_pdb, kabsch_rmsd, lddt_ca,
                            map_to_target, parse_tmscore_output, read_ca_structure,
                            run_tmscore, sha256_file)


def atom_line(serial, resname, number, xyz, chain="A", insertion=" ", atom="CA",
              altloc=" ", occupancy=1.0):
    x, y, z = xyz
    return (f"ATOM  {serial:5d} {atom:>4s}{altloc}{resname:3s} {chain}{number:4d}{insertion}   "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{10.0:6.2f}           C  ")


def fixture_structure(sequence, identifiers=None):
    identifiers = identifiers or [str(i + 1) for i in range(len(sequence))]
    residues = [{"residue_id": identifier, "amino_acid": aa,
                 "xyz": np.array([i, i * i, 0.0])}
                for i, (aa, identifier) in enumerate(zip(sequence, identifiers))]
    return {"sequence": sequence, "residues": residues}


def official_format_output(length=4, normalization=8, common=4, user_score="0.3750",
                           gdt="0.7500"):
    """Synthetic scores in official TMscore.h full-output format, not results."""
    return (f"Structure1: model.pdb    Length={length:5d}\n"
            f"Structure2: reference.pdb    Length={length:5d} (by which all scores are normalized)\n"
            f"Number of residues in common={common:5d}\n"
            "RMSD of  the common residues=    2.100\n\n"
            "TM-score    = 0.7000  (d0= 0.50)\n"
            f"GDT-TS-score= {gdt} %(d<1)=0.5000 %(d<2)=0.7500 %(d<4)=0.7500 %(d<8)=1.0000\n"
            f"TM-score    = {user_score}  (if normalized by user-specified LN={normalization:.2f} and d0=0.50)\n")


class RigidAndDistanceTests(unittest.TestCase):
    def setUp(self):
        # Noncoplanar, asymmetric geometry exposes improper reflection fits.
        self.reference = np.array([[0., 0., 0.], [3., 0., 0.], [0., 4., 0.], [0., 0., 5.]])

    def test_rigid_transform_invariance(self):
        angle = 0.61
        q = np.array([[np.cos(angle), -np.sin(angle), 0.],
                      [np.sin(angle), np.cos(angle), 0.], [0., 0., 1.]])
        mobile = self.reference @ q + [20., -3., 8.]
        fitted = kabsch_rmsd(mobile, self.reference)
        self.assertLess(fitted["rmsd_angstrom"], 1e-12)
        np.testing.assert_allclose(fitted["fitted"], self.reference, atol=1e-12)
        np.testing.assert_allclose(mobile @ fitted["rotation_row_vector"] + fitted["translation"],
                                   self.reference, atol=1e-12)
        self.assertAlmostEqual(lddt_ca(mobile, self.reference)["score"], 1.0)

    def test_perturbation_degrades_metrics_without_pruning(self):
        mobile = self.reference.copy()
        mobile[-1] += [4., 3., 2.]
        fitted = kabsch_rmsd(mobile, self.reference)
        self.assertGreater(fitted["rmsd_angstrom"], 0.5)
        self.assertEqual(len(fitted["distances_angstrom"]), 4)
        self.assertLess(lddt_ca(mobile, self.reference)["score"], 1.0)

    def test_reflection_is_not_accepted_as_rigid_rotation(self):
        mirrored = self.reference * [-1., 1., 1.]
        fitted = kabsch_rmsd(mirrored, self.reference)
        self.assertAlmostEqual(np.linalg.det(fitted["rotation_row_vector"]), 1.0)
        self.assertGreater(fitted["rmsd_angstrom"], 1.0)
        # Pairwise-distance scores alone do not detect chirality.
        self.assertAlmostEqual(lddt_ca(mirrored, self.reference)["score"], 1.0)

    def test_nonfinite_coordinates_rejected_for_fit(self):
        for value in (np.nan, np.inf, -np.inf):
            coordinates = self.reference.copy()
            coordinates[0, 0] = value
            with self.assertRaises(ValueError):
                kabsch_rmsd(coordinates, self.reference)

    def test_unpaired_or_too_few_coordinates_rejected(self):
        with self.assertRaises(ValueError):
            kabsch_rmsd(self.reference[:2], self.reference[:2])
        with self.assertRaises(ValueError):
            kabsch_rmsd(self.reference, self.reference[:3])

    def test_lddt_strict_threshold_and_reference_radius(self):
        reference = np.array([[0., 0., 0.], [1., 0., 0.]])
        mobile = np.array([[0., 0., 0.], [1.5, 0., 0.]])
        result = lddt_ca(mobile, reference)
        self.assertEqual(result["conserved_by_threshold"], [0, 1, 1, 1])
        self.assertEqual(result["score"], 0.75)
        # Inclusion is strictly '< radius', not '<='.
        self.assertIsNone(lddt_ca(reference, reference, radius=1.)["score"])

    def test_missing_mobile_row_penalizes_contacts(self):
        reference = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]])
        mobile = reference.copy()
        mobile[2] = np.nan
        result = lddt_ca(mobile, reference)
        self.assertEqual(result["eligible_reference_pairs"], 3)
        self.assertEqual(result["conserved_by_threshold"], [1, 1, 1, 1])
        self.assertAlmostEqual(result["score"], 1 / 3)
        self.assertEqual(result["missing_mobile_residues"], 1)

    def test_lddt_rejects_infinite_or_partial_missing_coordinates(self):
        for value in (np.inf, np.nan):
            mobile = self.reference.copy()
            mobile[0, 0] = value
            with self.assertRaises(ValueError):
                lddt_ca(mobile, self.reference)
        for radius in (0., -1., np.inf, np.nan):
            with self.assertRaises(ValueError):
                lddt_ca(self.reference, self.reference, radius)

    def test_no_reference_contacts_is_unavailable(self):
        far = np.array([[0., 0., 0.], [30., 0., 0.]])
        result = lddt_ca(far, far)
        self.assertIsNone(result["score"])
        self.assertEqual(result["status"], "NO_REFERENCE_CONTACTS")


class StructureAndMappingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / "coordinates.pdb"

    def write_lines(self, lines):
        self.path.write_text("\n".join(lines + ["END", ""]), encoding="ascii")

    def test_missing_ca_and_insertion_codes_preserve_sequence_anchor(self):
        self.write_lines([atom_line(1, "ALA", 10, [0., 0., 0.]),
                          atom_line(2, "ALA", 11, [1., 0., 0.], atom="N"),
                          atom_line(3, "CYS", 11, [2., 0., 0.], insertion="A"),
                          atom_line(4, "ALA", 12, [3., 0., 0.])])
        structure = read_ca_structure(self.path, "A")
        self.assertEqual(structure["sequence"], "ACA")
        self.assertEqual(structure["polymer_sequence"], "AACA")
        self.assertEqual(structure["missing_ca"][0]["residue_id"], "11")
        mapping = map_to_target("AACA", structure)
        self.assertEqual(list(mapping), [1, 3, 4])
        self.assertEqual(mapping[3]["residue_id"], "11A")
        self.assertEqual(structure["sha256"], sha256_file(self.path))

    def test_explicit_model_required_for_ensemble(self):
        self.write_lines(["MODEL        7", atom_line(1, "ALA", 1, [0., 0., 0.]),
                          "ENDMDL", "MODEL        9", atom_line(1, "ALA", 1, [5., 0., 0.]), "ENDMDL"])
        with self.assertRaisesRegex(MappingError, "select model_index"):
            read_ca_structure(self.path, "A")
        structure = read_ca_structure(self.path, "A", model_index=1)
        self.assertEqual(structure["model_serial_number"], 9)
        self.assertEqual(structure["residues"][0]["xyz"][0], 5.)
        for index in (-1, 2, 0.5, True):
            with self.assertRaises(MappingError):
                read_ca_structure(self.path, "A", model_index=index)

    def test_explicit_chain_selection_and_missing_chain(self):
        self.write_lines([atom_line(1, "ALA", 1, [0., 0., 0.]),
                          atom_line(2, "CYS", 1, [8., 0., 0.], chain="B")])
        self.assertEqual(read_ca_structure(self.path, "B")["sequence"], "C")
        for chain in ("Z", "", None):
            with self.assertRaises(MappingError):
                read_ca_structure(self.path, chain)

    def test_nonfinite_ca_rejected_at_read(self):
        self.write_lines([atom_line(1, "ALA", 1, [np.nan, 0., 0.])])
        with self.assertRaisesRegex(MappingError, "nonfinite"):
            read_ca_structure(self.path, "A")

    def test_tied_alternate_ca_requires_selection(self):
        self.write_lines([atom_line(1, "ALA", 1, [0., 0., 0.], altloc="A", occupancy=0.5),
                          atom_line(2, "ALA", 1, [9., 0., 0.], altloc="B", occupancy=0.5)])
        with self.assertRaisesRegex(MappingError, "tied alternate"):
            read_ca_structure(self.path, "A")
        structure = read_ca_structure(self.path, "A", altloc="B")
        self.assertEqual(structure["residues"][0]["altloc"], "B")
        self.assertEqual(structure["residues"][0]["xyz"][0], 9.)

    def test_mmcif_author_chain_is_used(self):
        pdb = "\n".join([atom_line(1, "ALA", 42, [0., 0., 0.]),
                          atom_line(2, "CYS", 42, [2., 0., 0.], insertion="A"), "END", ""])
        structure = PDBParser(QUIET=True).get_structure("test", StringIO(pdb))
        structure[0]["A"].id = "AUTHOR_LONG"
        cif = self.path.with_suffix(".cif")
        writer = MMCIFIO()
        writer.set_structure(structure)
        writer.save(str(cif))
        parsed = read_ca_structure(cif, "AUTHOR_LONG")
        self.assertEqual(parsed["chain"], "AUTHOR_LONG")
        self.assertEqual([r["residue_id"] for r in parsed["residues"]], ["42", "42A"])

    def test_ambiguous_sequence_alignment_is_blocked(self):
        with self.assertRaisesRegex(MappingError, "Ambiguous"):
            map_to_target("AAAA", fixture_structure("AAA"))
        mapped = map_to_target("AAAA", fixture_structure("AAA"), {2: "1", 3: "2", 4: "3"})
        self.assertEqual(list(mapped), [2, 3, 4])

    def test_sequence_mismatch_and_invalid_target_blocked(self):
        with self.assertRaisesRegex(MappingError, "mismatch"):
            map_to_target("ACDE", fixture_structure("ACNE"))
        for sequence in ("", "acde", "AC-X", "ACDX"):
            with self.assertRaises(MappingError):
                map_to_target(sequence, fixture_structure("ACDE"))

    def test_explicit_mapping_identity_unique_positions_and_residues(self):
        structure = fixture_structure("AAC", ["42", "42A", "43"])
        good = map_to_target("AAC", structure, {"1": "42", 2: "42A", 3: "43"})
        self.assertEqual(good[2]["residue_id"], "42A")
        invalid = [{1: "43"}, {1: "42", 2: "42"}, {1: "42", "1": "42A"},
                   {0: "42"}, {4: "42"}, {1.5: "42"}, {True: "42"}, {}, {1: "99"}]
        for mapping in invalid:
            with self.subTest(mapping=mapping), self.assertRaises(MappingError):
                map_to_target("AAC", structure, mapping)

    def test_fixed_mask_writer_roundtrip_and_validation(self):
        xyz = np.array([[1., 2., 3.], [4., 5., 6.], [7., 8., 9.]])
        _write_ca_pdb(self.path, [1, 3, 4], "ACDE", xyz)
        parsed = read_ca_structure(self.path, "A")
        self.assertEqual([r["number"] for r in parsed["residues"]], [1, 3, 4])
        self.assertEqual(parsed["sequence"], "ADE")
        np.testing.assert_allclose([r["xyz"] for r in parsed["residues"]], xyz)
        with self.assertRaises(ValueError):
            _write_ca_pdb(self.path, [1, 3], "ACDE", xyz)
        with self.assertRaises(ValueError):
            _write_ca_pdb(self.path, [1, 1, 4], "ACDE", xyz)
        with self.assertRaises(ValueError):
            _write_ca_pdb(self.path, [1, 3, 4], "ACDE", xyz * np.inf)
        with self.assertRaises(ValueError):
            _write_ca_pdb(self.path, [1, 3, 4], "ACDE", xyz * 1e6)


class OfficialTMscoreTests(unittest.TestCase):
    def test_real_official_executable_invariance_and_perturbation_when_available(self):
        executable = os.environ.get("AGENT2_TMSCORE_EXECUTABLE") or shutil.which("TMscore")
        if executable is None:
            self.skipTest("Official TMscore optional: set AGENT2_TMSCORE_EXECUTABLE or add to PATH")
        reference = np.array([[0., 0., 0.], [3., 0., 0.], [0., 4., 0.], [0., 0., 5.]])
        mobile = reference @ np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]]) + [10., 2., 7.]
        with tempfile.TemporaryDirectory() as directory:
            model_path, reference_path = Path(directory) / "model.pdb", Path(directory) / "ref.pdb"
            _write_ca_pdb(reference_path, [1, 2, 3, 4], "ACDE", reference)
            _write_ca_pdb(model_path, [1, 2, 3, 4], "ACDE", mobile)
            invariant = run_tmscore(model_path, reference_path, 8, 4, executable)
            self.assertEqual(invariant["status"], "AVAILABLE", invariant)
            self.assertAlmostEqual(invariant["tm_score"], 0.5, places=3)
            self.assertAlmostEqual(invariant["gdt_ts"], 0.5, places=3)
            mobile[-1] += [30., 0., 0.]
            _write_ca_pdb(model_path, [1, 2, 3, 4], "ACDE", mobile)
            perturbed = run_tmscore(model_path, reference_path, 8, 4, executable)
            self.assertEqual(perturbed["status"], "AVAILABLE", perturbed)
            self.assertLess(perturbed["tm_score"], invariant["tm_score"])
            self.assertLessEqual(perturbed["gdt_ts"], invariant["gdt_ts"])

    def test_full_output_selects_user_normalization_and_rescales_gdt(self):
        result = parse_tmscore_output(official_format_output(), reference_length=8, mask_count=4)
        self.assertEqual(result["tm_score"], 0.375)
        self.assertEqual(result["gdt_ts_native_mask_normalized"], 0.75)
        self.assertEqual(result["gdt_ts"], 0.375)
        self.assertEqual(result["reference_normalization_length"], 8)

    def test_wrong_normalization_or_correspondence_rejected(self):
        invalid = [official_format_output(normalization=9), official_format_output(common=3),
                   official_format_output(length=3), official_format_output(user_score="1.1000"),
                   official_format_output(user_score="0.7500"), official_format_output(gdt="-0.1")]
        for stdout in invalid:
            with self.subTest(stdout=stdout), self.assertRaises(ValueError):
                parse_tmscore_output(stdout, 8, 4)

    def test_tmalign_or_missing_user_normalization_rejected(self):
        tmalign = "TM-score= 0.50000 (normalized by length of Structure_2: L=8, d0=0.50)\n"
        with self.assertRaises(ValueError):
            parse_tmscore_output(tmalign, 8, 4)
        output = official_format_output().split("TM-score    = 0.3750")[0]
        with self.assertRaises(ValueError):
            parse_tmscore_output(output, 8, 4)
        with self.assertRaises(ValueError):
            parse_tmscore_output(official_format_output() * 2, 8, 4)

    def test_missing_gdt_is_null_not_zero(self):
        stdout = "\n".join(line for line in official_format_output().splitlines()
                           if not line.startswith("GDT-TS-score"))
        result = parse_tmscore_output(stdout, 8, 4)
        self.assertEqual(result["tm_score"], 0.375)
        self.assertIsNone(result["gdt_ts"])
        self.assertEqual(result["gdt_status"], "UNAVAILABLE_IN_OUTPUT")

    def test_invalid_normalization_inputs_rejected(self):
        for reference_length, mask_count in ((3, 4), (8, 2), (8.0, 4), (True, 4), (8, 4.0)):
            with self.assertRaises(ValueError):
                parse_tmscore_output(official_format_output(), reference_length, mask_count)

    def test_optional_executable_unavailable_does_not_fake_scores(self):
        with patch("agent2.metrics.shutil.which", return_value=None):
            result = run_tmscore("model.pdb", "reference.pdb", 8, 4,
                                 executable="definitely_not_an_installed_tmscore")
        self.assertEqual(result["status"], "UNAVAILABLE")
        self.assertIsNone(result["tm_score"])
        self.assertIsNone(result["gdt_ts"])

    def test_wrapper_uses_fixed_correspondence_and_keeps_provenance(self):
        completed = subprocess.CompletedProcess([], 0, official_format_output(), "")
        with patch("agent2.metrics.shutil.which", return_value="TMscore.exe"), \
                patch("agent2.metrics.subprocess.run", return_value=completed) as run, \
                patch("agent2.metrics.sha256_file", return_value="a" * 64):
            result = run_tmscore("model.pdb", "reference.pdb", 8, 4)
        self.assertEqual(result["tm_score"], 0.375)
        self.assertEqual(run.call_args.args[0][-2:], ["-l", "8"])
        self.assertNotIn("-seq", run.call_args.args[0])
        self.assertEqual(result["executable_sha256"], "a" * 64)

    def test_invalid_executable_output_keeps_null_and_raw_evidence(self):
        completed = subprocess.CompletedProcess([], 0, "not TMscore output", "diagnostic")
        with patch("agent2.metrics.shutil.which", return_value="TMscore.exe"), \
                patch("agent2.metrics.subprocess.run", return_value=completed):
            result = run_tmscore("model.pdb", "reference.pdb", 8, 4)
        self.assertEqual(result["status"], "UNAVAILABLE")
        self.assertIsNone(result["tm_score"])
        self.assertEqual(result["stdout"], "not TMscore output")
        self.assertEqual(result["stderr"], "diagnostic")


if __name__ == "__main__":
    unittest.main()
