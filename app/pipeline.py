"""Agent 1 demo pipeline for one target, run step by step into runs/<timestamp>/.

  python app/pipeline.py init    --run-dir runs/<ts> [--mode agent|direct]
  python app/pipeline.py step S  --run-dir runs/<ts>     S = search, select, align, model, compare
  python app/pipeline.py all     --run-dir runs/<ts>     init (if needed) and all five steps

Every step records its state, timing and a plain-language summary in
<run-dir>/status.json, appends to <run-dir>/pipeline.log, and writes tool
output to <run-dir>/logs/<step>.log. Saved results in targets/ are only read.
AlphaFold3 is not run: the model computed ahead of time in
targets/<T>/af3_model/ is compared with SCRIPTS/compare_to_af3.py.
"""

import argparse
import csv
import json
import math
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime

APP = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(APP)
SCRIPTS = os.path.join(REPO, "SCRIPTS")
sys.path.insert(0, SCRIPTS)

TARGET = "T1147"
START, END = 1, 103
PY = sys.executable
CPUS = int(os.environ.get("SLURM_CPUS_PER_TASK") or 4)

HH_CONDA = "/cluster/software/conda-envs-global/rosetta_1_2025/miniconda3"
UNIREF = "/cluster/VAST/cietest/rosettafold_databases/UniRef30_2020_06/UniRef30_2020_06"
PDB70 = "/cluster/VAST/cietest/alphafold_database/pdb70/pdb70"
AF3_DIR = os.path.join(REPO, "targets", TARGET, "af3_model", TARGET.lower())

STEPS = [
    ("search", "Search for templates"),
    ("select", "Choose a template"),
    ("align", "Align the sequence to the template"),
    ("model", "Build 3D models with MODELLER"),
    ("compare", "Compare with the AlphaFold3 model"),
]


# ---------------------------------------------------------------- status file

def status_path(run):
    return os.path.join(run, "status.json")


def load_status(run):
    with open(status_path(run)) as handle:
        return json.load(handle)


def save_status(run, status):
    tmp = status_path(run) + ".tmp"
    with open(tmp, "w") as handle:
        json.dump(status, handle, indent=2)
    os.replace(tmp, status_path(run))


def log(run, text):
    stamp = datetime.now().strftime("%H:%M:%S")
    with open(os.path.join(run, "pipeline.log"), "a") as handle:
        for line in text.rstrip().splitlines():
            handle.write(f"[{stamp}] {line}\n")


def init_run(run, mode):
    os.makedirs(os.path.join(run, "logs"), exist_ok=True)
    tdir = os.path.join(run, TARGET)
    os.makedirs(tdir, exist_ok=True)
    shutil.copy2(os.path.join(REPO, "targets", TARGET, f"{TARGET}.fasta"), os.path.join(tdir, f"{TARGET}.fasta"))
    status = {
        "target": TARGET, "mode": mode, "cpus": CPUS, "host": os.uname().nodename,
        "created": time.time(), "finished": None, "state": "pending",
        "steps": [{"id": sid, "title": title, "state": "pending", "start": None, "end": None,
                   "seconds": None, "summary": "", "error": ""} for sid, title in STEPS],
    }
    save_status(run, status)
    log(run, f"Run created for {TARGET} on {status['host']} with {CPUS} CPUs ({mode} mode).")


def run_cmd(run, step, cmd, shell=False, cwd=None):
    """Run a command, append its output to logs/<step>.log, raise on failure."""
    with open(os.path.join(run, "logs", f"{step}.log"), "a") as handle:
        handle.write(f"$ {cmd if shell else ' '.join(cmd)}\n")
        handle.flush()
        proc = subprocess.run(cmd, shell=shell, cwd=cwd, stdout=handle, stderr=subprocess.STDOUT,
                              executable="/bin/bash" if shell else None)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed (exit {proc.returncode}); see logs/{step}.log")


# ---------------------------------------------------------------- steps

def paths(run):
    t = os.path.join(run, TARGET)
    return {
        "dir": t,
        "fasta": os.path.join(t, f"{TARGET}.fasta"),
        "a3m": os.path.join(t, "hhsearch_output", f"{TARGET}.a3m"),
        "hhr": os.path.join(t, "hhsearch_output", f"{TARGET}.hhr"),
        "templates": os.path.join(t, "templates"),
        "pir": os.path.join(t, "alignment", f"{TARGET}.pir"),
        "models": os.path.join(t, "tbm_models"),
        "comparison": os.path.join(t, "comparison"),
        "viewer": os.path.join(t, "viewer"),
    }


def step_search(run, p):
    from top_hits import parse_hhr
    os.makedirs(os.path.dirname(p["a3m"]), exist_ok=True)
    activate = f"set +u; source {HH_CONDA}/etc/profile.d/conda.sh && conda activate RoseTTAFold; set -u"
    log(run, f"hhblits: 3 iterations against UniRef30 (June 2020), {CPUS} CPUs")
    run_cmd(run, "search", f"{activate}; hhblits -i {p['fasta']} -d {UNIREF} -oa3m {p['a3m']} "
                           f"-o /dev/null -n 3 -e 0.001 -cpu {CPUS}", shell=True)
    n_seqs = sum(1 for line in open(p["a3m"]) if line.startswith(">"))
    log(run, f"hhblits found {n_seqs} sequences. hhsearch: searching pdb70 (April 2020)")
    run_cmd(run, "search", f"{activate}; hhsearch -i {p['a3m']} -d {PDB70} -o {p['hhr']} -cpu {CPUS}",
            shell=True)
    _, hits = parse_hhr(p["hhr"])
    top = hits[0]
    return (f"Collected {n_seqs} related sequences and searched the PDB. Found {len(hits)} possible "
            f"templates; the strongest is {top['id']} ({top['prob']:.1f}% probability).")


def step_select(run, p):
    from top_hits import parse_hhr
    csv_out = os.path.join(run, "screening.csv")
    run_cmd(run, "select", [PY, os.path.join(SCRIPTS, "screen_templates.py"), TARGET,
                            "--targets-dir", run, "--csv", csv_out])
    row = next(csv.DictReader(open(csv_out)))
    if row["status"] != "PASS":
        raise RuntimeError(f"No template meets the cutoffs (top hit {row['top_hit_pdb']}, "
                           f"identity {row['identity']}%)")
    qlen, hits = parse_hhr(p["hhr"])
    hit = next(h for h in hits if h["id"] == row["best_template"])
    sel = {"template": hit["id"], "rank": hit["rank"], "probability": hit["prob"], "evalue": hit["evalue"],
           "identity": hit["identity"], "coverage": round(hit["q_coverage"], 1),
           "query_range": [hit["q_start"], hit["q_end"]], "protein": hit["description"].split(";")[0]}
    json.dump(sel, open(os.path.join(run, "selection.json"), "w"), indent=2)
    log(run, f"Selected {sel['template']} (hit {sel['rank']}): prob {sel['probability']}, "
             f"identity {sel['identity']}%, coverage {sel['coverage']}%")
    return (f"Picked {sel['template']} ({sel['protein']}). It passes all cutoffs: "
            f"{sel['probability']:.1f}% probability, {sel['identity']}% identical, "
            f"covers residues {sel['query_range'][0]}-{sel['query_range'][1]} ({sel['coverage']:.0f}%).")


def step_align(run, p):
    sel = json.load(open(os.path.join(run, "selection.json")))
    run_cmd(run, "align", [PY, os.path.join(SCRIPTS, "prepare_templates.py"), sel["template"],
                           "-o", p["templates"]])
    run_cmd(run, "align", [PY, os.path.join(SCRIPTS, "build_pir.py"), p["hhr"], p["fasta"],
                           "--hit", str(sel["rank"]), "-o", p["pir"]])
    text = open(os.path.join(run, "logs", "align.log")).read()
    line = [l for l in text.splitlines() if l.startswith("Hit ")][-1]
    pairs = line.split("aligned pairs ")[1].split()[0]
    ident = line.split("identical ")[1].split()[0]
    log(run, line)
    return (f"Downloaded and cleaned chain {sel['template']} and lined up the sequence with it: "
            f"{pairs} residue pairs, {ident} of them identical.")


def step_model(run, p):
    sel = json.load(open(os.path.join(run, "selection.json")))
    code = sel["template"].split("_")[0].lower() + "_" + sel["template"].split("_")[1]
    run_cmd(run, "model", [PY, os.path.join(SCRIPTS, "run_modeller.py"), "--alignment", p["pir"],
                           "--template", code, "--template-dir", p["templates"], "--query", TARGET,
                           "--start", str(START), "--end", str(END), "--n-models", "5",
                           "--out-dir", p["models"]])
    from finalize_target import tbm_scores
    scores = tbm_scores(os.path.join(run, "logs", "model.log"))
    json.dump([{"file": f, "dope": d, "ga341": g} for f, d, g in scores],
              open(os.path.join(run, "models.json"), "w"), indent=2)
    best = min(scores, key=lambda r: r[1])
    log(run, f"Best model {best[0]}: DOPE {best[1]:.1f}, GA341 {best[2]:.3f}")
    return (f"Built 5 models of residues {START}-{END}. The best by DOPE score is model "
            f"{int(best[0].split('.')[1][-1])} (DOPE {best[1]:.0f}, GA341 {best[2]:.2f}).")


def step_compare(run, p):
    from pymol import cmd
    from finalize_target import tmscore
    from compare_to_af3 import find_tmalign
    from visualize_model import dope_profile

    best = min(json.load(open(os.path.join(run, "models.json"))), key=lambda m: m["dope"])
    best_pdb = os.path.join(p["models"], best["file"])
    af3_model = os.path.join(AF3_DIR, f"{TARGET.lower()}_model.cif")
    log(run, "Loading the AlphaFold3 model computed ahead of time: "
             f"{os.path.relpath(af3_model, REPO)}")
    run_cmd(run, "compare", [PY, os.path.join(SCRIPTS, "compare_to_af3.py"), best_pdb, af3_model,
                             "--out-dir", p["comparison"]])
    pair_dir = os.path.join(p["comparison"], f"{os.path.splitext(best['file'])[0]}_vs_{TARGET.lower()}_model")
    pair = json.load(open(os.path.join(pair_dir, "summary.json")))

    # Residue-matched comparison over the modeled range
    os.makedirs(p["viewer"], exist_ok=True)
    af3_range = os.path.join(p["viewer"], "af3_range.pdb")
    cmd.reinitialize()
    cmd.load(af3_model, "af3")
    cmd.remove(f"not polymer.protein or hydro or not resi {START}-{END}")
    cmd.save(af3_range, "af3")
    matched = tmscore(find_tmalign(), best_pdb, af3_range)

    # Viewer files: TBM with DOPE in B, superposed TBM, AF3 with pLDDT in B
    dope = dope_profile(best_pdb, os.path.join(p["viewer"], "tbm_dope.profile"))
    cmd.reinitialize()
    cmd.load(best_pdb, "tbm")
    resv = []
    cmd.iterate("tbm and name CA", "resv_list.append(int(resv))", space={"resv_list": resv})
    dope_by = dict(zip(resv, dope))
    cmd.alter("tbm", "b = dope_by.get(int(resv), 0.0)", space={"dope_by": dope_by})
    cmd.save(os.path.join(p["viewer"], "tbm_dope.pdb"), "tbm")
    shutil.copy2(os.path.join(pair_dir, "tbm_superposed.pdb"), os.path.join(p["viewer"], "tbm_superposed.pdb"))
    shutil.copy2(os.path.join(pair_dir, "af3.pdb"), os.path.join(p["viewer"], "af3.pdb"))

    # Per-residue CA deviation by residue number under TM-align's superposition
    cmd.reinitialize()
    cmd.load(os.path.join(pair_dir, "tbm_superposed.pdb"), "t")
    cmd.load(os.path.join(pair_dir, "af3.pdb"), "a")
    t_xyz, a_xyz, plddt = {}, {}, {}
    cmd.iterate_state(1, "t and name CA", "t_xyz[int(resv)] = (x, y, z)", space={"t_xyz": t_xyz})
    cmd.iterate_state(1, "a and name CA", "a_xyz[int(resv)] = (x, y, z); plddt[int(resv)] = b",
                      space={"a_xyz": a_xyz, "plddt": plddt})
    with open(os.path.join(run, "per_residue.csv"), "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["resi", "ca_deviation_A", "tbm_dope", "af3_plddt"])
        for r in sorted(t_xyz):
            if r in a_xyz:
                writer.writerow([r, round(math.dist(t_xyz[r], a_xyz[r]), 3),
                                 round(dope_by.get(r, float("nan")), 4), round(plddt[r], 2)])

    conf = json.load(open(os.path.join(AF3_DIR, f"{TARGET.lower()}_summary_confidences.json")))
    sel = json.load(open(os.path.join(run, "selection.json")))
    results = {
        "target": TARGET, "template": sel["template"], "identity": sel["identity"],
        "coverage": sel["coverage"], "query_range": sel["query_range"], "modeled_range": [START, END],
        "best_model": os.path.relpath(best_pdb, run), "best_dope": best["dope"], "best_ga341": best["ga341"],
        "af3_model": os.path.relpath(af3_model, REPO), "af3_mean_plddt": pair["af3_mean_plddt_all"],
        "af3_ptm": conf.get("ptm"), "af3_ranking_score": conf.get("ranking_score"),
        "tm_align": {k: pair[k] for k in ("tm_score_norm_tbm", "tm_score_norm_af3", "rmsd",
                                          "aligned_length", "pairs_within_5A")},
        "residue_matched": matched,
    }
    json.dump(results, open(os.path.join(run, "results.json"), "w"), indent=2)
    log(run, f"TM-score {pair['tm_score_norm_tbm']:.3f} (TM-align), {matched['tm_score']:.3f} "
             f"(residue-matched); RMSD {pair['rmsd']:.2f} A over {pair['aligned_length']} residues")
    return (f"Compared the best model with the AlphaFold3 model computed ahead of time: TM-score "
            f"{pair['tm_score_norm_tbm']:.2f} (above 0.5 means the same fold), RMSD {pair['rmsd']:.1f} A "
            f"over {pair['aligned_length']} matched residues.")


STEP_FUNCS = {"search": step_search, "select": step_select, "align": step_align,
              "model": step_model, "compare": step_compare}


def run_step(run, sid):
    status = load_status(run)
    idx = [s["id"] for s in status["steps"]].index(sid)
    for prev in status["steps"][:idx]:
        if prev["state"] != "done":
            raise SystemExit(f"Step '{prev['id']}' has not finished; run the steps in order.")
    step = status["steps"][idx]
    step.update(state="running", start=time.time(), end=None, seconds=None, summary="", error="")
    status["state"] = "running"
    save_status(run, status)
    log(run, f"Step {idx + 1}/5: {step['title']}")
    try:
        summary = STEP_FUNCS[sid](run, paths(run))
    except Exception as exc:  # record the failure for the timeline, then exit non-zero
        status = load_status(run)
        step = status["steps"][idx]
        step.update(state="failed", end=time.time(), error=str(exc))
        step["seconds"] = round(step["end"] - step["start"], 1)
        status["state"] = "failed"
        status["finished"] = time.time()
        save_status(run, status)
        log(run, f"FAILED: {exc}")
        print(f"FAILED: {step['title']}: {exc}")
        raise SystemExit(1)
    status = load_status(run)
    step = status["steps"][idx]
    step.update(state="done", end=time.time(), summary=summary)
    step["seconds"] = round(step["end"] - step["start"], 1)
    if idx == len(STEPS) - 1:
        status["state"] = "done"
        status["finished"] = time.time()
    save_status(run, status)
    log(run, f"Done in {step['seconds']:.0f} s. {summary}")
    print(f"Step {idx + 1}/5 done in {step['seconds']:.0f} s: {summary}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("action", choices=["init", "step", "all"])
    parser.add_argument("step", nargs="?", choices=[s for s, _ in STEPS])
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--mode", default="direct", choices=["agent", "direct"])
    args = parser.parse_args()
    run = os.path.abspath(args.run_dir)

    if args.action == "init":
        init_run(run, args.mode)
    elif args.action == "step":
        if not args.step:
            parser.error("step needs a step name")
        run_step(run, args.step)
    else:
        if not os.path.exists(status_path(run)):
            init_run(run, args.mode)
        for sid, _ in STEPS:
            run_step(run, sid)


if __name__ == "__main__":
    main()
