"""Sequence-verified, fixed-correspondence C-alpha structure evaluation.

Adapted from this team's Project 1 coordinate-analysis helpers. RMSD uses an
independent NumPy Kabsch implementation. The lDDT formula follows Mariani et al.
(2013), https://doi.org/10.1093/bioinformatics/btt473, and the documented
OpenStructure distance tests; no OpenStructure implementation is copied.
TM/GDT are obtained only from an optional official TMscore executable. Its
USalign source has a permissive notice-based license, not a public-domain
designation: https://github.com/pylelab/USalign/blob/master/LICENSE.

Coordinates are in Angstroms. Callers must establish reference provenance and
freeze a common target-position mask before comparing multiple predictions.
Sequence matching alone does not prove that a reference is the CASP target.
"""
from __future__ import annotations

import hashlib
import math
from numbers import Integral
from pathlib import Path
import re
import shutil
import subprocess
from typing import Mapping

import numpy as np
from Bio import Align
from Bio.PDB import MMCIFParser, PDBParser, is_aa
from Bio.PDB.PDBExceptions import PDBConstructionException
from Bio.SeqUtils import seq1, seq3


AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
TMSCORE_SOURCE = "https://github.com/pylelab/USalign/blob/master/TMscore.cpp"
LDDT_SOURCE = "https://openstructure.org/docs/2.11/mol/alg/lddt/"
__all__ = ["MappingError", "read_ca_structure", "map_to_target", "kabsch_rmsd",
           "lddt_ca", "parse_tmscore_output", "run_tmscore", "_write_ca_pdb", "sha256_file"]


class MappingError(ValueError):
    """Selection or sequence correspondence cannot be established safely."""


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


_sha = sha256_file


def _select_ca(residue, altloc: str | None):
    atom = residue["CA"]
    if atom.is_disordered() != 2:
        return atom
    alternatives = list(atom.child_dict.values())
    if altloc is not None:
        if altloc not in atom.child_dict:
            raise MappingError(f"Residue {residue.id}: CA alternate location {altloc!r} absent")
        return atom.child_dict[altloc]
    occupancies = [a.get_occupancy() for a in alternatives]
    if any(x is None or not math.isfinite(float(x)) for x in occupancies):
        raise MappingError(f"Residue {residue.id}: alternate CA occupancy unresolved; select altloc")
    maximum = max(occupancies)
    best = [a for a in alternatives if a.get_occupancy() == maximum]
    if len(best) != 1:
        raise MappingError(f"Residue {residue.id}: tied alternate CA occupancies; select altloc")
    return best[0]


def read_ca_structure(path: str | Path, chain: str, model_index: int | None = None,
                      *, altloc: str | None = None) -> dict:
    """Read one explicit author chain/model from PDB or mmCIF.

    model_index is zero-based, distinct from a PDB MODEL serial number. Multiple
    models require explicit selection. Insertion codes and author numbering are
    retained. Missing CA residues are recorded and retain their sequence anchors.
    Alternate CA atoms use a unique highest occupancy; ties require altloc.
    Ambiguous residue identities and nonfinite selected coordinates are rejected.
    """
    path = Path(path)
    if not isinstance(chain, str) or not chain:
        raise MappingError("An explicit, nonempty author chain ID is required")
    if altloc is not None and (not isinstance(altloc, str) or len(altloc) != 1):
        raise MappingError("altloc must be a single-character alternate-location ID")
    if model_index is not None and (isinstance(model_index, bool) or not isinstance(model_index, Integral)):
        raise MappingError("model_index must be a zero-based integer")
    if path.suffix.lower() in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True, auth_chains=True, auth_residues=True)
    else:
        parser = PDBParser(QUIET=True, PERMISSIVE=False)
    try:
        structure = parser.get_structure("coordinates", str(path))
    except (ValueError, KeyError, PDBConstructionException) as exc:
        raise MappingError(f"{path.name}: invalid coordinate file: {exc}") from exc
    models = list(structure.get_models())
    if not models:
        raise MappingError(f"{path.name}: no coordinate models")
    if len(models) != 1 and model_index is None:
        raise MappingError(f"{path.name}: {len(models)} coordinate models; select model_index explicitly")
    index = 0 if model_index is None else int(model_index)
    if not 0 <= index < len(models):
        raise MappingError(f"Invalid model index {index} for {path.name}")
    if chain not in models[index]:
        raise MappingError(f"{path.name}: chain {chain!r} absent; present: {[c.id for c in models[index]]}")
    residues, missing_ca, polymer_residues = [], [], []
    seen = set()
    for residue in models[index][chain]:
        if residue.is_disordered() == 2:
            raise MappingError(f"{path.name}: alternative residue identities at {residue.id}; clean explicitly")
        if not is_aa(residue, standard=False):
            continue
        number, insertion = int(residue.id[1]), residue.id[2].strip()
        identifier = f"{number}{insertion}"
        if identifier in seen:
            raise MappingError("Duplicate author residue IDs require an explicitly cleaned coordinate copy")
        seen.add(identifier)
        amino_acid = seq1(residue.resname, custom_map={"MSE": "M", "SEC": "U", "PYL": "O"})
        entry = {"chain": chain, "number": number, "insertion_code": insertion,
                 "residue_id": identifier, "resname": residue.resname, "amino_acid": amino_acid}
        polymer_residues.append(entry)
        if "CA" not in residue:
            missing_ca.append(entry)
            continue
        atom = _select_ca(residue, altloc)
        xyz = np.asarray(atom.coord, dtype=float)
        if xyz.shape != (3,) or not np.isfinite(xyz).all():
            raise MappingError(f"{path.name}: nonfinite CA coordinates at author residue {identifier}")
        entry.update({"xyz": xyz.copy(), "ca_bfactor": float(atom.bfactor),
                      "altloc": atom.get_altloc().strip(), "occupancy": atom.get_occupancy()})
        residues.append(entry)
    if not residues:
        raise MappingError(f"{path.name}: selected chain has no protein C-alpha coordinates")
    return {"path": str(path.resolve()), "sha256": sha256_file(path), "chain": chain,
            "model_index": index, "model_serial_number": getattr(models[index], "serial_num", None),
            "residues": residues, "missing_ca": missing_ca,
            "sequence": "".join(r["amino_acid"] for r in residues),
            "polymer_residues": polymer_residues,
            "polymer_sequence": "".join(r["amino_acid"] for r in polymer_residues),
            "altloc_policy": "explicit " + altloc if altloc is not None else "unique highest occupancy"}


def _target_sequence(target: str) -> None:
    if not isinstance(target, str) or not target or set(target) - AA:
        raise MappingError("Target sequence must be nonempty uppercase canonical protein without gaps")


def _position(value) -> int:
    if isinstance(value, bool) or not (isinstance(value, Integral) or
                                      isinstance(value, str) and re.fullmatch(r"[1-9][0-9]*", value)):
        raise MappingError(f"Target position {value!r} must be a positive integer")
    return int(value)


def map_to_target(target: str, structure: dict,
                  explicit_mapping: Mapping | None = None) -> dict[int, dict]:
    """Map 1-based target positions to sequence-identical residues with CA atoms.

    Explicit values are author residue IDs, e.g. {1: '42', 2: '42A'}. Residue
    reuse, duplicate target positions (including 1 and '1'), and identity
    mismatches are rejected. Positions omitted from the mapping stay unscored.
    Automatic mapping requires a unique global sequence alignment and never
    substitutes a structure alignment. Known missing-CA residues anchor sequence
    alignment, but are excluded from the returned coordinate correspondence.
    """
    _target_sequence(target)
    residues = structure["residues"]
    if not residues:
        raise MappingError("Structure has no observed CA residues")
    ids = [str(r["residue_id"]) for r in residues]
    if len(ids) != len(set(ids)):
        raise MappingError("Structure contains duplicate author residue IDs")
    if explicit_mapping is not None:
        if not isinstance(explicit_mapping, Mapping):
            raise MappingError("Explicit residue mapping must be a target-position dictionary")
        lookup = dict(zip(ids, residues))
        result, used = {}, set()
        for raw_position, raw_identifier in explicit_mapping.items():
            position, identifier = _position(raw_position), str(raw_identifier)
            if (not 1 <= position <= len(target) or identifier not in lookup or
                    identifier in used or position in result):
                raise MappingError(f"Invalid/repeated explicit mapping: target {position} -> residue {identifier}")
            residue = lookup[identifier]
            if residue["amino_acid"] != target[position - 1]:
                raise MappingError(f"Explicit mapping sequence mismatch at target position {position}")
            result[position] = residue
            used.add(identifier)
        if not result:
            raise MappingError("Explicit residue mapping is empty")
        return dict(sorted(result.items()))
    polymer = structure.get("polymer_residues", residues)
    observed = "".join(r["amino_acid"] for r in polymer)
    if not observed:
        raise MappingError("Structure has no protein sequence")
    if observed == target:
        pairs = [(i, i) for i in range(len(target))]
    else:
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score, aligner.mismatch_score = 2.0, -3.0
        aligner.open_gap_score, aligner.extend_gap_score = -5.0, -0.5
        alignments = iter(aligner.align(target, observed))
        first = next(alignments, None)
        if first is None:
            raise MappingError("No sequence alignment")
        if next(alignments, None) is not None:
            raise MappingError("Ambiguous optimal sequence-to-coordinate mapping; provide explicit mapping")
        pairs = [(ti, ri) for (t0, t1), (r0, r1) in zip(first.aligned[0], first.aligned[1])
                 for ti, ri in zip(range(int(t0), int(t1)), range(int(r0), int(r1)))]
    result = {}
    for ti, ri in pairs:
        if target[ti] != observed[ri]:
            raise MappingError(f"Sequence mismatch at target {ti + 1}; construct differences require explicit mapping")
        residue = polymer[ri]
        if "xyz" in residue:
            result[ti + 1] = residue
    if not result:
        raise MappingError("No sequence-identical CA residue correspondences")
    return result


def _paired_arrays(mobile, reference, minimum: int = 0):
    mobile, reference = np.asarray(mobile, dtype=float), np.asarray(reference, dtype=float)
    if (mobile.shape != reference.shape or mobile.ndim != 2 or mobile.shape[1] != 3
            or len(mobile) < minimum):
        raise ValueError(f"Paired Nx3 arrays with at least {minimum} points are required")
    return mobile, reference


def kabsch_rmsd(mobile: np.ndarray, reference: np.ndarray) -> dict:
    """Least-squares rigid fit on every supplied pair; determinant +1, no pruning."""
    mobile, reference = _paired_arrays(mobile, reference, 3)
    if not np.isfinite(mobile).all() or not np.isfinite(reference).all():
        raise ValueError("Coordinates must be finite")
    m_center, r_center = mobile.mean(axis=0), reference.mean(axis=0)
    u, _, vt = np.linalg.svd((mobile - m_center).T @ (reference - r_center))
    correction = np.eye(3)
    correction[-1, -1] = -1.0 if np.linalg.det(u @ vt) < 0 else 1.0
    rotation = u @ correction @ vt
    translation = r_center - m_center @ rotation
    fitted = mobile @ rotation + translation
    distances = np.linalg.norm(fitted - reference, axis=1)
    return {"rmsd_angstrom": float(np.sqrt(np.mean(distances ** 2))),
            "rotation_row_vector": rotation, "translation": translation,
            "fitted": fitted, "distances_angstrom": distances}


def lddt_ca(mobile: np.ndarray, reference: np.ndarray, radius: float = 15.0) -> dict:
    """Superposition-free C-alpha lDDT on a caller-supplied fixed residue mask.

    Unordered distinct-residue reference pairs with distance < radius are tested
    with strict errors < 0.5, 1, 2, 4 A. This is not all-atom lDDT and performs no
    stereochemical checks. All-NaN mobile rows denote missing coordinates and
    fail their reference contacts; infinities and partly missing rows are invalid.
    An empty reference contact set returns None, never a fabricated zero.
    """
    mobile, reference = _paired_arrays(mobile, reference)
    if not np.isfinite(reference).all() or not math.isfinite(radius) or radius <= 0:
        raise ValueError("Reference coordinates must be finite and radius finite and positive")
    complete = np.isfinite(mobile).all(axis=1)
    missing = np.isnan(mobile).all(axis=1)
    if not np.all(complete | missing):
        raise ValueError("Mobile coordinates must be finite or entirely NaN per missing residue")
    ref_dist = np.linalg.norm(reference[:, None] - reference[None, :], axis=2)
    mob_dist = np.linalg.norm(mobile[:, None] - mobile[None, :], axis=2)
    eligible = np.triu(ref_dist < radius, k=1)
    errors = np.abs(mob_dist[eligible] - ref_dist[eligible])
    conserved = [int(np.count_nonzero(errors < t)) for t in (0.5, 1.0, 2.0, 4.0)]
    pair_count = int(eligible.sum())
    return {"score": sum(conserved) / (4 * pair_count) if pair_count else None,
            "eligible_reference_pairs": pair_count, "conserved_by_threshold": conserved,
            "thresholds_angstrom": [0.5, 1.0, 2.0, 4.0], "radius_angstrom": float(radius),
            "missing_mobile_residues": int(missing.sum()),
            "scope": "lDDT-Calpha on supplied fixed residue mask; no stereochemical checks",
            "status": "AVAILABLE" if pair_count else "NO_REFERENCE_CONTACTS", "source": LDDT_SOURCE}


def _normalization(reference_length: int, mask_count: int) -> None:
    if (isinstance(reference_length, bool) or not isinstance(reference_length, Integral) or
            isinstance(mask_count, bool) or not isinstance(mask_count, Integral) or
            not reference_length >= mask_count >= 3):
        raise ValueError("Reference normalization length must be integer >= common mask size >= 3")


def parse_tmscore_output(stdout: str, reference_length: int, mask_count: int) -> dict:
    """Parse official full TMscore output with explicit -l normalization.

    The official optimized GDT-TS output is divided by Structure2's coordinate
    count even when -l is used for TM-score; rescale it to the declared reference
    length. Both input lengths and the common count must equal the frozen mask.
    TMalign output and guessed normalization are rejected.
    """
    _normalization(reference_length, mask_count)
    numeric = r"([0-9.eE+-]+)"
    scores = []
    for raw, explanation in re.findall(r"^\s*TM-score\s*=\s*" + numeric + r"\s*\(([^\n]+)\)", stdout, re.MULTILINE):
        if "normalized by user-specified" not in explanation:
            continue
        length = re.search(r"\bLN?\s*=\s*" + numeric, explanation)
        if not length or not math.isclose(float(length.group(1)), reference_length, rel_tol=0, abs_tol=1e-8):
            raise ValueError("TMscore output normalization does not match declared reference length")
        scores.append(float(raw))
    if len(scores) != 1:
        raise ValueError("TMscore output must contain one user-specified reference-length score (-l)")
    for name in ("Structure1", "Structure2"):
        lengths = re.findall(r"^\s*" + name + r":.*?Length=\s*(\d+)", stdout, re.MULTILINE)
        if len(lengths) != 1 or int(lengths[0]) != mask_count:
            raise ValueError("TMscore input length does not match the frozen residue correspondence")
    common = re.findall(r"Number of residues in common=\s*(\d+)", stdout)
    if len(common) != 1 or int(common[0]) != mask_count:
        raise ValueError("TMscore did not preserve the entire frozen residue correspondence")
    gdts = re.findall(r"^\s*GDT-TS-score\s*=\s*" + numeric, stdout, re.MULTILINE)
    if len(gdts) > 1:
        raise ValueError("TMscore output contains multiple GDT-TS scores")
    gdt_native = float(gdts[0]) if gdts else None
    tm_score = scores[0]
    if (not math.isfinite(tm_score) or not 0 <= tm_score <= 1 or
            gdt_native is not None and (not math.isfinite(gdt_native) or not 0 <= gdt_native <= 1)):
        raise ValueError("TMscore returned a score outside [0,1]")
    # The published executable prints four decimals, hence a half-unit tolerance.
    if tm_score > mask_count / reference_length + 0.000051:
        raise ValueError("TM-score exceeds the maximum possible for this mask and normalization")
    return {"tm_score": tm_score,
            "gdt_ts": gdt_native * mask_count / reference_length if gdt_native is not None else None,
            "gdt_ts_native_mask_normalized": gdt_native,
            "reference_normalization_length": int(reference_length),
            "frozen_mask_count": int(mask_count),
            "gdt_normalization": "Official optimized GDT counts / declared reference length (rescaled from native mask length)",
            "gdt_status": "AVAILABLE" if gdt_native is not None else "UNAVAILABLE_IN_OUTPUT",
            "status": "AVAILABLE", "source": TMSCORE_SOURCE}


def run_tmscore(model_path: str | Path, reference_path: str | Path, reference_length: int,
                mask_count: int, executable: str | Path | None = None) -> dict:
    """Run optional official TMscore on already-renumbered fixed-mask CA PDBs.

    No executable download, shell interpolation, sequence alignment, or proxy
    score. Missing tools or failed/invalid output return null TM/GDT with evidence.
    Input/normalization programming errors raise before executable discovery.
    """
    _normalization(reference_length, mask_count)
    requested = str(executable) if executable is not None else "TMscore"
    resolved = shutil.which(requested)
    if resolved is None and Path(requested).is_file():
        resolved = str(Path(requested).resolve())
    requirement = "Provide official pylelab/USalign TMscore (not TMalign); see " + TMSCORE_SOURCE
    unavailable = {"status": "UNAVAILABLE", "tm_score": None, "gdt_ts": None,
                   "reference_normalization_length": int(reference_length), "requirement": requirement}
    if not resolved:
        return {**unavailable, "diagnostic": f"Executable not found: {requested}"}
    command = [resolved, str(Path(model_path).resolve()), str(Path(reference_path).resolve()),
               "-l", str(reference_length)]
    completed = None
    try:
        completed = subprocess.run(command, capture_output=True, text=True, timeout=120, check=True)
        parsed = parse_tmscore_output(completed.stdout, reference_length, mask_count)
        return {**parsed, "command": command, "executable_sha256": sha256_file(resolved),
                "stdout": completed.stdout, "stderr": completed.stderr}
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        return {**unavailable, "diagnostic": str(exc), "command": command,
                "stdout": getattr(completed, "stdout", None) or getattr(exc, "stdout", None),
                "stderr": getattr(completed, "stderr", None) or getattr(exc, "stderr", None)}


def _write_ca_pdb(path: str | Path, positions: list[int], sequence: str, xyz: np.ndarray) -> None:
    """Write the frozen mask with target residue numbers shared by both inputs.

    Insertion codes are deliberately absent after correspondence has been frozen.
    This PDB is a scoring artifact, not a replacement for the original coordinate
    file. PDB limits are checked so columns never overflow silently.
    """
    _target_sequence(sequence)
    normalized = [_position(p) for p in positions]
    coords = np.asarray(xyz, dtype=float)
    if not normalized or len(normalized) > 99999 or coords.shape != (len(normalized), 3):
        raise ValueError("A nonempty coordinate array matching the mask is required")
    if normalized != sorted(set(normalized)):
        raise ValueError("Mask positions must be unique and strictly increasing")
    if not all(1 <= p <= min(len(sequence), 9999) for p in normalized):
        raise ValueError("Mask positions exceed the target sequence or PDB residue-number limit")
    if not np.isfinite(coords).all():
        raise ValueError("Coordinates must be finite")
    lines = []
    for serial, (position, point) in enumerate(zip(normalized, coords), start=1):
        fields = [f"{value:8.3f}" for value in point]
        if any(len(field) != 8 for field in fields):
            raise ValueError("Coordinates exceed the PDB field width")
        residue = seq3(sequence[position - 1]).upper()
        lines.append(f"ATOM  {serial:5d}  CA  {residue:3s} A{position:4d}    "
                     + "".join(fields) + f"{1.0:6.2f}{0.0:6.2f}           C  ")
    Path(path).write_text("\n".join(lines + ["TER", "END", ""]), encoding="ascii")
