"""Compare a template-based model with an AlphaFold3 model using TM-align.

Accepts PDB or mmCIF files (AlphaFold Server and local AlphaFold3 write
mmCIF). Both structures are reduced to protein atoms and written as PDB,
aligned with TM-align, and the TBM model is superposed onto the AF3 model
with TM-align's rotation matrix. Outputs go to <out-dir>/<tbm>_vs_<af3>/:

  tmalign_output.txt         raw TM-align report
  tmalign_matrix.txt         TM-align rotation/translation (TBM -> AF3)
  summary.json, summary.txt  TM-scores, RMSD, aligned length, AF3 pLDDT
  per_residue_distances.csv  CA-CA distance for each aligned residue pair
  superposition.png          TBM (blue) on AF3 (orange; unaligned AF3 grey)
  superposition.pse          PyMOL session of the same scene
"""

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess

from pymol import cmd

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
TMALIGN_FALLBACK = "/cluster/pixstor/dkhf3-lab/users/nancy/envs/tbm_tools/bin/TMalign"


def find_tmalign(path=None):
    for candidate in (path, shutil.which("TMalign"), TMALIGN_FALLBACK):
        if candidate and os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
    raise FileNotFoundError("TMalign not found; activate tbm_tools or pass --tmalign")


def to_clean_pdb(src, dst, name):
    """Write the protein atoms of src (PDB or mmCIF) to dst as PDB."""
    cmd.load(src, name)
    cmd.remove(f"{name} and not polymer.protein")
    cmd.remove(f"{name} and hydro")
    cmd.save(dst, name)
    cmd.delete(name)


def ca_records(name):
    """Return [(chain, resi, resn, (x, y, z), b)] for CA atoms in file order."""
    records = []
    cmd.iterate_state(1, f"{name} and name CA and not alt B",
                      "records.append((chain, resi, resn, (x, y, z), b))",
                      space={"records": records})
    return records


def run_tmalign(tmalign, pdb1, pdb2, matrix_file):
    result = subprocess.run([tmalign, pdb1, pdb2, "-m", matrix_file],
                            capture_output=True, text=True, check=True)
    return result.stdout


def parse_tmalign(text):
    """Extract scores and the three-line alignment from TM-align stdout."""
    m = re.search(r"Aligned length=\s*(\d+), RMSD=\s*([\d.]+), Seq_ID=n_identical/n_aligned=\s*([\d.]+)", text)
    # Older TM-align builds say "Chain_N", newer ones "Structure_N"
    tms = dict((c, float(s)) for s, c in re.findall(
        r"TM-score=\s*([\d.]+) \((?:if )?normalized by length of (?:Chain|Structure)_(\d)", text))
    lengths = dict(re.findall(r"Length of (?:Chain|Structure)_(\d):\s*(\d+) residues", text))
    lines = text.rstrip().splitlines()
    marker = next(i for i, line in enumerate(lines) if line.startswith('(":" denotes'))
    seq1, match, seq2 = lines[marker + 1], lines[marker + 2], lines[marker + 3]
    return {
        "aligned_length": int(m.group(1)),
        "rmsd": float(m.group(2)),
        "seq_identity_aligned": float(m.group(3)),
        "tm_score_norm_tbm": tms["1"],
        "tm_score_norm_af3": tms["2"],
        "length_tbm": int(lengths["1"]),
        "length_af3": int(lengths["2"]),
    }, (seq1, match, seq2)


def read_matrix(path):
    """Return (t, U) from a TM-align -m file; x' = t + U x."""
    rows = []
    with open(path) as handle:
        for line in handle:
            parts = line.split()
            if len(parts) == 5 and parts[0] in ("0", "1", "2"):
                rows.append([float(v) for v in parts[1:]])
    t = [r[0] for r in rows]
    U = [r[1:] for r in rows]
    return t, U


def aligned_pairs(alignment, ca1, ca2):
    """Map alignment columns to (index1, index2, marker) residue index pairs."""
    seq1, match, seq2 = alignment
    pairs, i, j = [], 0, 0
    for a, mk, b in zip(seq1, match, seq2):
        if a != "-" and b != "-":
            pairs.append((i, j, mk))
        i += a != "-"
        j += b != "-"
    if i != len(ca1) or j != len(ca2):
        raise ValueError("TM-align alignment does not match CA atoms read by PyMOL")
    return pairs


def transform(xyz, t, U):
    return tuple(t[k] + sum(U[k][n] * xyz[n] for n in range(3)) for k in range(3))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("tbm", help="Template-based model (PDB or mmCIF)")
    parser.add_argument("af3", help="AlphaFold3 model (PDB or mmCIF)")
    parser.add_argument("--out-dir", default="comparison_results")
    parser.add_argument("--tmalign", default=None, help="Path to TMalign binary")
    parser.add_argument("--width", type=int, default=1600)
    parser.add_argument("--height", type=int, default=1200)
    args = parser.parse_args()

    tmalign = find_tmalign(args.tmalign)
    tbm_stem = os.path.splitext(os.path.basename(args.tbm))[0]
    af3_stem = os.path.splitext(os.path.basename(args.af3))[0]
    out = os.path.join(args.out_dir, f"{tbm_stem}_vs_{af3_stem}")
    os.makedirs(out, exist_ok=True)

    # Clean inputs to protein-only PDB
    cmd.reinitialize()
    tbm_pdb = os.path.join(out, "tbm.pdb")
    af3_pdb = os.path.join(out, "af3.pdb")
    to_clean_pdb(os.path.abspath(args.tbm), tbm_pdb, "tbm_in")
    to_clean_pdb(os.path.abspath(args.af3), af3_pdb, "af3_in")

    # TM-align: structure 1 = TBM, structure 2 = AF3
    matrix_file = os.path.join(out, "tmalign_matrix.txt")
    report = run_tmalign(tmalign, tbm_pdb, af3_pdb, matrix_file)
    with open(os.path.join(out, "tmalign_output.txt"), "w") as handle:
        handle.write(report)
    scores, alignment = parse_tmalign(report)
    t, U = read_matrix(matrix_file)

    # Load cleaned structures and superpose TBM onto AF3 with TM-align's matrix
    cmd.load(tbm_pdb, "tbm")
    cmd.load(af3_pdb, "af3")
    ttt = [U[0][0], U[0][1], U[0][2], t[0],
           U[1][0], U[1][1], U[1][2], t[1],
           U[2][0], U[2][1], U[2][2], t[2],
           0.0, 0.0, 0.0, 1.0]
    cmd.transform_selection("tbm", ttt, homogenous=1)
    cmd.save(os.path.join(out, "tbm_superposed.pdb"), "tbm")

    ca_tbm, ca_af3 = ca_records("tbm"), ca_records("af3")
    pairs = aligned_pairs(alignment, ca_tbm, ca_af3)

    # Per-residue CA distances after superposition
    dist_rows = []
    for i, j, mk in pairs:
        d = math.dist(ca_tbm[i][3], ca_af3[j][3])
        dist_rows.append({
            "tbm_chain": ca_tbm[i][0], "tbm_resi": ca_tbm[i][1], "tbm_resn": ca_tbm[i][2],
            "af3_chain": ca_af3[j][0], "af3_resi": ca_af3[j][1], "af3_resn": ca_af3[j][2],
            "ca_distance": round(d, 3), "tmalign_marker": mk,
        })
    with open(os.path.join(out, "per_residue_distances.csv"), "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dist_rows[0]))
        writer.writeheader()
        writer.writerows(dist_rows)

    # Sanity check: RMSD over aligned pairs under TM-align's superposition
    rmsd_check = math.sqrt(sum(r["ca_distance"] ** 2 for r in dist_rows) / len(dist_rows))

    # AF3 writes pLDDT into the B-factor column
    plddt_all = [r[4] for r in ca_af3]
    plddt_aligned = [ca_af3[j][4] for _, j, _ in pairs]
    scores.update({
        "tbm_file": os.path.abspath(args.tbm),
        "af3_file": os.path.abspath(args.af3),
        "rmsd_aligned_pairs_tm_superposition": round(rmsd_check, 3),
        "af3_mean_plddt_all": round(sum(plddt_all) / len(plddt_all), 2),
        "af3_mean_plddt_aligned": round(sum(plddt_aligned) / len(plddt_aligned), 2),
        "pairs_within_5A": sum(mk == ":" for _, _, mk in pairs),
    })
    with open(os.path.join(out, "summary.json"), "w") as handle:
        json.dump(scores, handle, indent=2)

    # Superposition image: TBM blue, aligned AF3 orange, unaligned AF3 grey
    aligned_af3 = {(ca_af3[j][0], ca_af3[j][1]) for _, j, _ in pairs}
    cmd.bg_color("white")
    cmd.hide("everything")
    cmd.show("cartoon")
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("ray_shadows", 0)
    cmd.set("antialias", 2)
    cmd.set_color("tbm_blue", [0x2a / 255, 0x78 / 255, 0xd6 / 255])
    cmd.set_color("af3_orange", [0xeb / 255, 0x68 / 255, 0x34 / 255])
    cmd.color("grey80", "af3")
    for chain, resi in aligned_af3:
        cmd.color("af3_orange", f"af3 and chain {chain} and resi {resi}")
    cmd.color("tbm_blue", "tbm")
    cmd.orient("tbm or af3")
    cmd.zoom("tbm or af3", 2)
    cmd.set("ray_opaque_background", 1)
    png = os.path.join(out, "superposition.png")
    cmd.png(png, width=args.width, height=args.height, dpi=300, ray=1)
    cmd.save(os.path.join(out, "superposition.pse"))

    lines = [
        f"TBM model : {scores['tbm_file']} ({scores['length_tbm']} residues)",
        f"AF3 model : {scores['af3_file']} ({scores['length_af3']} residues)",
        f"TM-score (normalized by TBM length) : {scores['tm_score_norm_tbm']:.4f}",
        f"TM-score (normalized by AF3 length) : {scores['tm_score_norm_af3']:.4f}",
        f"RMSD (TM-align, aligned residues)   : {scores['rmsd']:.2f} A",
        f"Aligned length                      : {scores['aligned_length']}",
        f"Aligned pairs within 5 A            : {scores['pairs_within_5A']}",
        f"Sequence identity of aligned pairs  : {scores['seq_identity_aligned']:.3f}",
        f"AF3 mean pLDDT (all / aligned)      : {scores['af3_mean_plddt_all']:.1f} / "
        f"{scores['af3_mean_plddt_aligned']:.1f}",
        f"Results: {out}",
    ]
    with open(os.path.join(out, "summary.txt"), "w") as handle:
        handle.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
