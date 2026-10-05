"""Finish one target after its MODELLER and AlphaFold3 jobs: compare and record.

Picks the best TBM model by DOPE from targets/<T>/tbm_models.log, runs
compare_to_af3.py against the top-ranked AF3 model, adds a residue-matched
TM-score (TMscore) over the modeled range, and writes targets/<T>/notes.txt
and targets/<T>/metadata.json.
"""

import argparse
import datetime
import json
import os
import re
import subprocess
import sys

from pymol import cmd

from compare_to_af3 import find_tmalign
from top_hits import parse_hhr

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
TBM = os.path.dirname(SCRIPTS)


def tbm_scores(log_file):
    """Return [(model_file, dope, ga341)] from the run_modeller.py summary."""
    rows, in_summary = [], False
    with open(log_file) as handle:
        for line in handle:
            if line.startswith("=== Summary ==="):
                in_summary = True
            elif in_summary:
                m = re.match(r"^(\S+\.pdb)\s+(-?[\d.]+)\s+([\d.]+)\s*$", line)
                if m:
                    rows.append((m.group(1), float(m.group(2)), float(m.group(3))))
    if not rows:
        raise ValueError(f"No model scores found in {log_file}")
    return rows


def tmscore(tmalign, model, reference):
    exe = os.path.join(os.path.dirname(tmalign), "TMscore")
    out = subprocess.run([exe, model, reference], capture_output=True, text=True, check=True).stdout
    grab = lambda pat: float(re.search(pat, out).group(1))
    return {"tm_score": grab(r"TM-score\s*=\s*([\d.]+)"),
            "rmsd": grab(r"RMSD of\s+the common residues=\s*([\d.]+)"),
            "gdt_ts": grab(r"GDT-TS-score=\s*([\d.]+)"),
            "residues_in_common": int(grab(r"Number of residues in common=\s*(\d+)"))}


def save_range(src, dst, start, end):
    cmd.delete("all")
    cmd.load(src, "obj")
    cmd.remove(f"not polymer.protein or hydro or not resi {start}-{end}")
    cmd.save(dst, "obj")
    cmd.delete("all")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target")
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    parser.add_argument("--hits", type=int, nargs="+", required=True, help="hhsearch hit numbers used")
    parser.add_argument("--template-choice", default="", help="Why these templates were chosen")
    parser.add_argument("--note", action="append", default=[], help="Problem or remark for notes.txt")
    parser.add_argument("--modeller-job", default="")
    parser.add_argument("--af3-job", default="")
    args = parser.parse_args()

    t = args.target
    d = os.path.join(TBM, "targets", t)
    tmalign = find_tmalign()

    # Inputs and templates
    fasta = os.path.join(d, f"{t}.fasta")
    header = open(fasta).readline()[1:].strip()
    length = len("".join(l.strip() for l in open(fasta) if not l.startswith(">")))
    qlen, hits = parse_hhr(os.path.join(d, "hhsearch_output", f"{t}.hhr"))
    by_rank = {h["rank"]: h for h in hits}
    templates = []
    for rank in args.hits:
        h = by_rank[rank]
        res = re.search(r"([\d.]+)A \{", h["description"])
        templates.append({
            "pdb_chain": h["id"], "hhsearch_rank": rank, "probability": h["prob"], "evalue": h["evalue"],
            "identity_pct": h["identity"], "query_range": [h["q_start"], h["q_end"]],
            "template_seqres_range": [h["t_start"], h["t_end"]],
            "resolution_A": float(res.group(1)) if res else None,
            "protein": h["description"].split(";")[0],
            "file": f"templates/{h['id'].split('_')[0].lower()}_{h['id'].split('_')[1]}.pdb"})

    # TBM models
    scores = tbm_scores(os.path.join(d, "tbm_models.log"))
    best_file, best_dope, best_ga341 = min(scores, key=lambda r: r[1])
    best_path = os.path.join(d, "tbm_models", best_file)

    # AF3 top model
    af3_dir = os.path.join(d, "af3_model", t.lower())
    af3_model = os.path.join(af3_dir, f"{t.lower()}_model.cif")
    conf = json.load(open(os.path.join(af3_dir, f"{t.lower()}_summary_confidences.json")))

    # Comparison (sequence-independent, TM-align) via compare_to_af3.py
    comp_dir = os.path.join(d, "comparison")
    subprocess.run([sys.executable, os.path.join(SCRIPTS, "compare_to_af3.py"), best_path, af3_model,
                    "--out-dir", comp_dir], check=True, stdout=subprocess.DEVNULL)
    pair_dir = os.path.join(comp_dir, f"{os.path.splitext(best_file)[0]}_vs_{t.lower()}_model")
    pair = json.load(open(os.path.join(pair_dir, "summary.json")))

    # Residue-matched comparison over the modeled range
    af3_range = os.path.join(pair_dir, f"af3_{args.start}-{args.end}.pdb")
    save_range(af3_model, af3_range, args.start, args.end)
    matched = tmscore(tmalign, best_path, af3_range)
    cmd.load(af3_model, "af3")
    plddt = {}
    cmd.iterate("af3 and name CA", "plddt[int(resv)] = b", space={"plddt": plddt})
    cmd.delete("all")
    in_range = [v for r, v in plddt.items() if args.start <= r <= args.end]

    problems = list(args.note)
    uncovered = [r for r in range(args.start, args.end + 1)
                 if not any(tp["query_range"][0] <= r <= tp["query_range"][1] for tp in templates)]
    if uncovered:
        spans, s0 = [], uncovered[0]
        for a, b in zip(uncovered, uncovered[1:] + [None]):
            if b != a + 1:
                spans.append(f"{s0}-{a}" if a != s0 else f"{s0}")
                s0 = b
        noun = "Residue" if len(uncovered) == 1 else "Residues"
        verb = "has" if len(uncovered) == 1 else "have"
        problems.append(f"{noun} {', '.join(spans)} {verb} no template; "
                        f"MODELLER built {'it' if len(uncovered) == 1 else 'them'} without template restraints.")
    if conf.get("has_clash"):
        problems.append("AF3 top model flagged has_clash.")

    meta = {
        "target": t, "description": header, "length": length, "fasta": f"{t}.fasta",
        "created": datetime.date.today().isoformat(),
        "template_search": {"hhblits_db": "UniRef30_2020_06", "hhsearch_db": "pdb70 (2020-04)",
                            "hhr": f"hhsearch_output/{t}.hhr", "choice": args.template_choice},
        "templates": templates,
        "modeled_range": [args.start, args.end],
        "alignment": f"alignment/{t}.pir",
        "modeller": {"version": "10.8", "n_models": len(scores), "slurm_job": args.modeller_job,
                     "models": [{"file": f"tbm_models/{f}", "dope": dp, "ga341": ga} for f, dp, ga in scores],
                     "best_model": f"tbm_models/{best_file}", "best_dope": best_dope, "best_ga341": best_ga341},
        "alphafold3": {"module": "alphafold3/alphafold3_v301_deepmind", "slurm_job": args.af3_job,
                       "model": os.path.relpath(af3_model, d), "seeds": 1, "samples": 5,
                       "ptm": conf.get("ptm"), "ranking_score": conf.get("ranking_score"),
                       "fraction_disordered": conf.get("fraction_disordered"),
                       "mean_plddt": pair["af3_mean_plddt_all"],
                       "mean_plddt_modeled_range": round(sum(in_range) / len(in_range), 2)},
        "comparison": {
            "dir": os.path.relpath(pair_dir, d),
            "tm_align": {"tm_score_norm_tbm": pair["tm_score_norm_tbm"],
                         "tm_score_norm_af3": pair["tm_score_norm_af3"], "rmsd": pair["rmsd"],
                         "aligned_length": pair["aligned_length"], "pairs_within_5A": pair["pairs_within_5A"]},
            "residue_matched": matched},
        "problems": problems,
    }
    with open(os.path.join(d, "metadata.json"), "w") as handle:
        json.dump(meta, handle, indent=2)

    lines = [f"{t}: {header}", f"Length {length}; modeled residues {args.start}-{args.end}", "",
             "Templates:"]
    for tp in templates:
        lines.append(f"  {tp['pdb_chain']} (hhsearch hit {tp['hhsearch_rank']}): prob {tp['probability']:.1f}, "
                     f"E {tp['evalue']:.2g}, identity {tp['identity_pct']}%, query {tp['query_range'][0]}-"
                     f"{tp['query_range'][1]}; {tp['protein']}"
                     + (f", {tp['resolution_A']} A" if tp["resolution_A"] else ""))
    if args.template_choice:
        lines.append(f"  Choice: {args.template_choice}")
    lines += ["", f"Best TBM model: tbm_models/{best_file} (DOPE {best_dope:.1f}, GA341 {best_ga341:.3f})",
              f"AF3 top model: {os.path.relpath(af3_model, d)} (mean pLDDT {pair['af3_mean_plddt_all']:.1f}, "
              f"pTM {conf.get('ptm')})",
              f"TBM vs AF3: TM-score {matched['tm_score']:.3f} residue-matched, "
              f"{pair['tm_score_norm_tbm']:.3f} TM-align; RMSD {matched['rmsd']:.2f} A residue-matched, "
              f"{pair['rmsd']:.2f} A TM-align over {pair['aligned_length']} residues",
              "", "Problems:"] + ([f"  - {p}" for p in problems] or ["  none"])
    with open(os.path.join(d, "notes.txt"), "w") as handle:
        handle.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
