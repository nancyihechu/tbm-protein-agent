"""Quality comparison of the TBM (MODELLER) models against AlphaFold3.

Runs compare_to_af3.py for the best TBM model vs the top-ranked AF3 model,
then scores every TBM model and every AF3 sample on the same residue range
with the same metrics, and writes presentation figures and tables to
<out-dir>/presentation/:

  fig1_superposition.png        TBM model on AF3 model (from compare_to_af3.py)
  fig2_af3_plddt.png            AF3 model colored by pLDDT
  fig3_per_residue.png          AF3 pLDDT, DOPE profiles, CA deviation per residue
  fig4_model_scores.png         z-DOPE and TM-score for every model
  model_scores.csv              one row per model (TBM 1-5, AF3 samples)
  per_residue.csv               per-residue pLDDT, DOPE and CA deviation
  summary.json, summary.md      headline numbers for the slides

DOPE and z-DOPE are computed for the AF3 model trimmed to the TBM range so
both methods are judged on the same residues.
"""

import argparse
import csv
import glob
import json
import math
import os
import re
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from modeller import Environ, Selection, log
from modeller.scripts import complete_pdb
from pymol import cmd

from compare_to_af3 import find_tmalign, parse_tmalign, run_tmalign

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
EXAMPLE = os.path.join(os.path.dirname(SCRIPTS), "examples", "test_target")

# Reference palette: categorical slots 1-2, text and surface tokens
TBM_COLOR = "#2a78d6"
AF3_COLOR = "#eb6834"
TBM_LIGHT = "#86b6ef"
TEXT = "#0b0b0b"
TEXT_2 = "#52514e"
GRID = "#e4e3df"
SURFACE = "#fcfcfb"
REGION = "#f0efec"
# AlphaFold pLDDT convention
PLDDT_BANDS = [(90, 100, "#0053d6", "Very high (>90)"), (70, 90, "#65cbf3", "Confident (70-90)"),
               (50, 70, "#ffdb13", "Low (50-70)"), (0, 50, "#ff7d45", "Very low (<50)")]


def env_for_dope():
    log.none()
    env = Environ()
    env.libs.topology.read(file="$(LIB)/top_heav.lib")
    env.libs.parameters.read(file="$(LIB)/par.lib")
    return env


def dope_scores(env, pdb, profile=None, window=15):
    """Return (DOPE, z-DOPE, per-residue profile or None) for a PDB file."""
    mdl = complete_pdb(env, pdb)
    sel = Selection(mdl)
    if profile:
        dope = sel.assess_dope(output="ENERGY_PROFILE NO_REPORT", file=profile,
                               normalize_profile=True, smoothing_window=window)
    else:
        dope = sel.assess_dope(output="NO_REPORT")
    zdope = mdl.assess_normalized_dope()
    values = None
    if profile:
        with open(profile) as handle:
            values = [float(l.split()[-1]) for l in handle if l.strip() and not l.startswith("#")]
    return dope, zdope, values


def save_protein_pdb(src, dst, resi_range=None):
    """Write protein atoms of src (PDB/mmCIF) to dst, optionally limited to a residue range."""
    cmd.delete("all")
    cmd.load(src, "obj")
    cmd.remove("not polymer.protein or hydro")
    if resi_range:
        cmd.remove(f"not resi {resi_range[0]}-{resi_range[1]}")
    cmd.save(dst, "obj")
    cmd.delete("obj")


def ca_values(pdb, prop="b"):
    """Return {resi: value} for CA atoms of a structure file."""
    cmd.delete("all")
    cmd.load(pdb, "obj")
    out = {}
    cmd.iterate_state(1, "obj and name CA and not alt B", f"out[int(resv)] = {prop}", space={"out": out})
    cmd.delete("obj")
    return out


def tm_pair(tmalign, pdb1, pdb2, workdir, tag):
    report = run_tmalign(tmalign, pdb1, pdb2, os.path.join(workdir, f"{tag}_matrix.txt"))
    scores, _ = parse_tmalign(report)
    return scores


def tmscore_seq(tmalign, model, reference):
    """Residue-matched TM-score/RMSD/GDT-TS (TMscore program), normalized by the reference."""
    tmscore = os.path.join(os.path.dirname(tmalign), "TMscore")
    out = subprocess.run([tmscore, model, reference], capture_output=True, text=True, check=True).stdout
    grab = lambda pat: float(re.search(pat, out).group(1))
    return {"tm": grab(r"TM-score\s*=\s*([\d.]+)"),
            "rmsd": grab(r"RMSD of\s+the common residues=\s*([\d.]+)"),
            "gdt_ts": grab(r"GDT-TS-score=\s*([\d.]+)"),
            "common": int(grab(r"Number of residues in common=\s*(\d+)"))}


def style_axes(ax):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(TEXT_2)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=TEXT_2, labelsize=9)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def render_plddt(af3_model, png, width=1600, height=1200):
    cmd.reinitialize()
    cmd.load(af3_model, "af3")
    cmd.remove("not polymer.protein or hydro")
    for lo, hi, hexcol, _ in PLDDT_BANDS:
        name = f"plddt_{lo}"
        cmd.set_color(name, [int(hexcol[i:i + 2], 16) / 255 for i in (1, 3, 5)])
        cmd.color(name, f"af3 and b > {lo - 0.001} and b < {hi + 0.001}")
    cmd.bg_color("white")
    cmd.hide("everything")
    cmd.show("cartoon")
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("ray_shadows", 0)
    cmd.set("antialias", 2)
    cmd.orient("af3")
    cmd.zoom("af3", 2)
    cmd.set("ray_opaque_background", 1)
    cmd.png(png, width=width, height=height, dpi=300, ray=1)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--af3-dir", default=os.path.join(EXAMPLE, "af3/output/test_target"))
    parser.add_argument("--tbm-dir", default=os.path.join(EXAMPLE, "models/first_pass"))
    parser.add_argument("--tbm-best", default="test_target.B99990003.pdb")
    parser.add_argument("--template", default=os.path.join(EXAMPLE, "templates/4ehx_A.pdb"))
    parser.add_argument("--start", type=int, default=103)
    parser.add_argument("--end", type=int, default=460)
    parser.add_argument("--regions", default="182-215,224-230,331-345")
    parser.add_argument("--window", type=int, default=15)
    parser.add_argument("--out-dir", default=os.path.join(EXAMPLE, "comparison_results"))
    args = parser.parse_args()

    pres = os.path.join(args.out_dir, "presentation")
    work = os.path.join(pres, "work")
    os.makedirs(work, exist_ok=True)
    tmalign = find_tmalign()
    regions = [tuple(int(v) for v in r.split("-")) for r in args.regions.split(",")]
    rng = (args.start, args.end)

    name = os.path.basename(os.path.normpath(args.af3_dir))
    af3_top = os.path.join(args.af3_dir, f"{name}_model.cif")
    tbm_best = os.path.join(args.tbm_dir, args.tbm_best)

    # 1. Full pairwise comparison of the headline pair
    subprocess.run([sys.executable, os.path.join(SCRIPTS, "compare_to_af3.py"), tbm_best, af3_top,
                    "--out-dir", args.out_dir], check=True, stdout=subprocess.DEVNULL)
    pair_dir = os.path.join(args.out_dir, f"{os.path.splitext(args.tbm_best)[0]}_{'vs'}_{name}_model")
    with open(os.path.join(pair_dir, "summary.json")) as handle:
        pair = json.load(handle)

    # 2. Collect models: TBM 1-5 and AF3 samples (+ top) as protein-only PDBs
    env = env_for_dope()
    with open(os.path.join(args.af3_dir, f"{name}_summary_confidences.json")) as handle:
        top_conf = json.load(handle)
    ranking = {}
    with open(os.path.join(args.af3_dir, "ranking_scores.csv")) as handle:
        for row in csv.DictReader(handle):
            ranking[f"seed-{row['seed']}_sample-{row['sample']}"] = float(row["ranking_score"])

    af3_top_full = os.path.join(work, "af3_top_full.pdb")
    af3_top_trim = os.path.join(work, "af3_top_trim.pdb")
    save_protein_pdb(af3_top, af3_top_full)
    save_protein_pdb(af3_top, af3_top_trim, rng)

    rows = []
    tbm_files = sorted(glob.glob(os.path.join(args.tbm_dir, "*.B9999*.pdb")))
    for pdb in tbm_files:
        label = "TBM " + os.path.basename(pdb).split(".")[1][-1]
        dope, zdope, _ = dope_scores(env, pdb)
        tm = tm_pair(tmalign, pdb, af3_top_full, work, label.replace(" ", "_"))
        tt = tm_pair(tmalign, pdb, args.template, work, label.replace(" ", "_") + "_tmpl")
        seq = tmscore_seq(tmalign, pdb, af3_top_trim)
        rows.append({"model": label, "method": "TBM (MODELLER)", "file": pdb,
                     "tmscore_seq_vs_af3_top": seq["tm"], "rmsd_seq_vs_af3_top": seq["rmsd"],
                     "gdt_ts_vs_af3_top": seq["gdt_ts"],
                     "residues": f"{args.start}-{args.end}", "dope": round(dope, 1),
                     "zdope": round(zdope, 3), "tm_vs_af3_top": tm["tm_score_norm_tbm"],
                     "rmsd_vs_af3_top": tm["rmsd"], "tm_vs_template": tt["tm_score_norm_tbm"],
                     "mean_plddt": "", "ptm": "", "ranking_score": "",
                     "best_tbm": os.path.basename(pdb) == args.tbm_best})

    for sample_dir in sorted(glob.glob(os.path.join(args.af3_dir, "seed-*_sample-*"))):
        tag = os.path.basename(sample_dir)
        cif = os.path.join(sample_dir, "model.cif")
        trim = os.path.join(work, f"af3_{tag}_trim.pdb")
        full = os.path.join(work, f"af3_{tag}_full.pdb")
        save_protein_pdb(cif, full)
        save_protein_pdb(cif, trim, rng)
        dope, zdope, _ = dope_scores(env, trim)
        plddt = ca_values(full)
        with open(os.path.join(sample_dir, "summary_confidences.json")) as handle:
            conf = json.load(handle)
        tm = tm_pair(tmalign, tbm_best, full, work, f"best_vs_{tag}")
        tt = tm_pair(tmalign, trim, args.template, work, f"{tag}_tmpl")
        rows.append({"model": f"AF3 {tag}", "method": "AlphaFold3", "file": cif,
                     "residues": f"{args.start}-{args.end} (trimmed)", "dope": round(dope, 1),
                     "zdope": round(zdope, 3), "tm_vs_af3_top": "", "rmsd_vs_af3_top": "",
                     "tm_vs_template": tt["tm_score_norm_tbm"],
                     "tm_best_tbm_vs_this": tm["tm_score_norm_tbm"],
                     "mean_plddt": round(sum(plddt.values()) / len(plddt), 2),
                     "ptm": conf.get("ptm"), "ranking_score": round(ranking.get(tag, float("nan")), 4),
                     "best_tbm": False})

    fields = ["model", "method", "residues", "dope", "zdope", "tmscore_seq_vs_af3_top",
              "rmsd_seq_vs_af3_top", "gdt_ts_vs_af3_top", "tm_vs_af3_top", "rmsd_vs_af3_top",
              "tm_best_tbm_vs_this", "tm_vs_template", "mean_plddt", "ptm", "ranking_score", "file"]
    with open(os.path.join(pres, "model_scores.csv"), "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    # 3. Per-residue data for the headline pair
    _, tbm_z, tbm_prof = dope_scores(env, tbm_best, os.path.join(work, "tbm_best.profile"), args.window)
    af3_dope, af3_z, af3_prof = dope_scores(env, af3_top_trim, os.path.join(work, "af3_top_trim.profile"),
                                            args.window)
    tbm_dope = next(r["dope"] for r in rows if r["best_tbm"])
    tbm_resi = sorted(ca_values(tbm_best))
    af3_trim_resi = sorted(ca_values(af3_top_trim))
    tbm_dope_by = dict(zip(tbm_resi, tbm_prof))
    af3_dope_by = dict(zip(af3_trim_resi, af3_prof))
    plddt_by = ca_values(af3_top_full)
    # Same sequence, so compare residue-by-residue (by number) under TM-align's superposition
    tbm_xyz = ca_values(os.path.join(pair_dir, "tbm_superposed.pdb"), "(x, y, z)")
    af3_xyz = ca_values(af3_top_full, "(x, y, z)")
    dev_by = {r: math.dist(tbm_xyz[r], af3_xyz[r]) for r in tbm_xyz if r in af3_xyz}
    with open(os.path.join(pair_dir, "per_residue_distances.csv")) as handle:
        aln = list(csv.DictReader(handle))
    in_register = sum(row["tbm_resi"] == row["af3_resi"] for row in aln)
    all_resi = sorted(plddt_by)
    with open(os.path.join(pres, "per_residue.csv"), "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["resi", "af3_plddt", "tbm_dope_smoothed", "af3_dope_smoothed", "ca_deviation_A"])
        for r in all_resi:
            writer.writerow([r, plddt_by.get(r, ""), tbm_dope_by.get(r, ""), af3_dope_by.get(r, ""),
                             dev_by.get(r, "")])

    def region_stats(r0, r1):
        sel = [r for r in range(r0, r1 + 1)]
        mean = lambda d: (sum(d[r] for r in sel if r in d) / max(1, sum(r in d for r in sel)))
        return {"region": f"{r0}-{r1}", "af3_plddt": round(mean(plddt_by), 1),
                "tbm_dope": round(mean(tbm_dope_by), 4), "af3_dope": round(mean(af3_dope_by), 4),
                "ca_deviation": round(mean(dev_by), 2),
                "resolved_pairs": sum(r in dev_by for r in sel), "length": r1 - r0 + 1}

    region_rows = [region_stats(*r) for r in [(1, args.start - 1)] + regions]

    # 4. Figures
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "text.color": TEXT,
                         "axes.labelcolor": TEXT, "figure.facecolor": SURFACE})

    render_plddt(af3_top_full, os.path.join(pres, "fig2_af3_plddt.png"))
    # Legend strip for the pLDDT render
    from matplotlib.patches import Patch
    fig = plt.figure(figsize=(7.5, 0.45))
    fig.legend(handles=[Patch(color=col, label=lab) for _, _, col, lab in PLDDT_BANDS],
               loc="center", ncol=4, frameon=False, fontsize=9, title="AlphaFold3 pLDDT",
               title_fontsize=9, handlelength=1.4, columnspacing=1.6)
    fig.savefig(os.path.join(pres, "fig2_plddt_legend.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    # fig3: three stacked panels, shared residue axis, one y-scale each
    fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True,
                             gridspec_kw={"height_ratios": [1, 1.2, 1], "hspace": 0.18})
    for ax in axes:
        style_axes(ax)
        ax.axvspan(1, args.start - 0.5, color=REGION, zorder=0)
        for r0, r1 in regions:
            ax.axvspan(r0 - 0.5, r1 + 0.5, color=REGION, zorder=0)
    axes[0].text((1 + args.start) / 2, 0.06, "no template\n(not in TBM model)",
                 transform=axes[0].get_xaxis_transform(), ha="center", va="bottom", fontsize=9, color=TEXT_2)
    for k, (r0, r1) in enumerate(regions):
        axes[0].text((r0 + r1) / 2, 0.06 + 0.13 * (k % 2), f"{r0}-{r1}",
                     transform=axes[0].get_xaxis_transform(), ha="center", va="bottom",
                     fontsize=9, color=TEXT_2)

    ax = axes[0]
    for lo, hi, col, _ in PLDDT_BANDS:
        ax.axhspan(lo, hi, color=col, alpha=0.10, zorder=0)
    ax.plot(all_resi, [plddt_by[r] for r in all_resi], color=AF3_COLOR, linewidth=2)
    ax.set_ylim(40, 100)
    ax.set_ylabel("AF3 pLDDT")
    ax.set_title("AlphaFold3 per-residue confidence", loc="left", fontsize=11, color=TEXT)

    ax = axes[1]
    ax.plot(tbm_resi, tbm_prof, color=TBM_COLOR, linewidth=2, label="TBM (MODELLER)")
    ax.plot(af3_trim_resi, af3_prof, color=AF3_COLOR, linewidth=2, label="AlphaFold3")
    ax.set_ylabel("DOPE per residue\n(smoothed, lower = better)")
    ax.set_title(f"Per-residue DOPE, residues {args.start}-{args.end}", loc="left", fontsize=11, color=TEXT)
    ax.legend(loc="lower left", frameon=False, fontsize=9, ncol=2)
    ax.text(tbm_resi[-1] + 3, tbm_prof[-1], "TBM", color=TEXT, fontsize=9, va="center")
    ax.text(af3_trim_resi[-1] + 3, af3_prof[-1], "AF3", color=TEXT, fontsize=9, va="center")

    ax = axes[2]
    dev_resi = sorted(dev_by)
    ax.plot(dev_resi, [dev_by[r] for r in dev_resi], color=TEXT_2, linewidth=2)
    ax.axhline(5, color=TEXT_2, linewidth=0.8, linestyle=(0, (4, 3)))
    ax.text(args.start, 5.3, "5 Å", fontsize=9, color=TEXT_2)
    ax.set_ylabel("CA deviation (Å)")
    ax.set_xlabel("Residue")
    ax.set_title("TBM vs AF3 deviation, same residue number (TM-align superposition)", loc="left", fontsize=11, color=TEXT)
    ax.set_xlim(1, all_resi[-1] + 12)
    fig.savefig(os.path.join(pres, "fig3_per_residue.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    # fig4: per-model z-DOPE (dot plot) and TM-score of each TBM model vs AF3
    tbm_rows = [r for r in rows if r["method"].startswith("TBM")]
    af3_rows = [r for r in rows if r["method"] == "AlphaFold3"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"wspace": 0.35})
    for ax in (a1, a2):
        style_axes(ax)
    labels = [r["model"] for r in tbm_rows] + [r["model"].replace("AF3 seed-", "AF3 s").replace("_sample-", "/")
                                                for r in af3_rows]
    vals = [r["zdope"] for r in tbm_rows + af3_rows]
    cols = [TBM_COLOR] * len(tbm_rows) + [AF3_COLOR] * len(af3_rows)
    y = list(range(len(vals)))[::-1]
    a1.hlines(y, [0] * len(vals), vals, color=GRID, linewidth=2, zorder=1)
    a1.scatter(vals, y, s=64, c=cols, edgecolors=SURFACE, linewidths=2, zorder=3)
    for yi, v in zip(y, vals):
        a1.text(v + 0.09, yi, f"{v:.2f}", va="center", ha="left", fontsize=8, color=TEXT)
    a1.set_yticks(y)
    a1.set_yticklabels(labels, fontsize=9)
    a1.axvline(0, color=TEXT_2, linewidth=0.8)
    a1.axvline(-1, color=TEXT_2, linewidth=0.8, linestyle=(0, (4, 3)))
    a1.text(-1.04, -0.95, "native-like < -1", fontsize=8, color=TEXT_2, ha="right", va="center")
    a1.set_ylim(-1.4, len(vals) - 0.4)
    a1.set_xlim(min(vals + [-1.2]) - 0.3, max(vals + [0.2]) + 0.45)
    a1.set_xlabel(f"z-DOPE, residues {args.start}-{args.end} (lower = better)")
    a1.xaxis.grid(True, color=GRID, linewidth=0.8)
    a1.yaxis.grid(False)
    a1.set_title("Model quality (normalized DOPE)", loc="left", fontsize=11)

    xs = list(range(len(tbm_rows)))
    series = [("Residue-matched (TMscore)", "tmscore_seq_vs_af3_top", TBM_COLOR, -0.19),
              ("Sequence-independent (TM-align)", "tm_vs_af3_top", TBM_LIGHT, 0.19)]
    for lab, key, col, dx in series:
        vals2 = [r[key] for r in tbm_rows]
        a2.bar([x + dx for x in xs], vals2, width=0.36, color=col, edgecolor=SURFACE,
               linewidth=2, label=lab)
        for x, v in zip(xs, vals2):
            a2.text(x + dx, v + 0.01, f"{v:.2f}", ha="center", fontsize=7, color=TEXT)
    a2.legend(loc="upper left", frameon=False, fontsize=8, bbox_to_anchor=(0, 1.0))
    a2.axhline(0.5, color=TEXT_2, linewidth=0.8, linestyle=(0, (4, 3)))
    a2.text(len(xs) - 0.4, 0.5, "same fold\n> 0.5", fontsize=8, color=TEXT_2, ha="left", va="center")
    a2.set_xlim(-0.6, len(xs) + 0.45)
    a2.set_xticks(xs)
    a2.set_xticklabels([r["model"] for r in tbm_rows], fontsize=9)
    a2.set_ylim(0, 1.1)
    a2.set_ylabel(f"TM-score vs AF3 top model\n(residues {args.start}-{args.end})")
    a2.set_title("Agreement with AlphaFold3", loc="left", fontsize=11)
    fig.savefig(os.path.join(pres, "fig4_model_scores.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    # fig1: copy the superposition render
    import shutil
    shutil.copy2(os.path.join(pair_dir, "superposition.png"), os.path.join(pres, "fig1_superposition.png"))

    # 5. Summary
    best_af3 = max(af3_rows, key=lambda r: r["ranking_score"])
    summary = {
        "tbm_model": tbm_best, "af3_model": af3_top,
        "tm_score_norm_tbm": pair["tm_score_norm_tbm"], "tm_score_norm_af3": pair["tm_score_norm_af3"],
        "rmsd": pair["rmsd"], "aligned_length": pair["aligned_length"],
        "pairs_within_5A": pair["pairs_within_5A"], "aligned_pairs_in_register": in_register,
        "tmscore_residue_matched": next(r["tmscore_seq_vs_af3_top"] for r in rows if r["best_tbm"]),
        "rmsd_residue_matched_tmscore": next(r["rmsd_seq_vs_af3_top"] for r in rows if r["best_tbm"]),
        "gdt_ts_residue_matched": next(r["gdt_ts_vs_af3_top"] for r in rows if r["best_tbm"]),
        "residue_matched_rmsd": round(math.sqrt(sum(d * d for d in dev_by.values()) / len(dev_by)), 2),
        "length_tbm": pair["length_tbm"],
        "length_af3": pair["length_af3"],
        "af3_ptm": top_conf.get("ptm"), "af3_ranking_score": top_conf.get("ranking_score"),
        "af3_fraction_disordered": top_conf.get("fraction_disordered"),
        "af3_has_clash": top_conf.get("has_clash"),
        "af3_mean_plddt_all": pair["af3_mean_plddt_all"], "af3_mean_plddt_aligned": pair["af3_mean_plddt_aligned"],
        "dope_tbm": tbm_dope, "dope_af3_trimmed": round(af3_dope, 1),
        "zdope_tbm": round(tbm_z, 3), "zdope_af3_trimmed": round(af3_z, 3),
        "tm_tbm_vs_template": next(r["tm_vs_template"] for r in rows if r["best_tbm"]),
        "tm_af3_vs_template": best_af3["tm_vs_template"],
        "regions": region_rows,
    }
    with open(os.path.join(pres, "summary.json"), "w") as handle:
        json.dump(summary, handle, indent=2)

    md = [
        "# TBM vs AlphaFold3: model quality comparison",
        "",
        f"- TBM model: `{os.path.relpath(tbm_best, SCRIPTS)}` (MODELLER, template 4EHX_A, residues {args.start}-{args.end})",
        f"- AF3 model: `{os.path.relpath(af3_top, SCRIPTS)}` (local AlphaFold3 3.0.0, full length 1-{pair['length_af3']}, top-ranked of {len(af3_rows)} samples)",
        "",
        f"## Structural agreement, residue-matched (TMscore, residues {args.start}-{args.end})",
        "",
        "Same sequence, so each TBM residue is compared only with the same AF3 residue.",
        "",
        "| Metric | Value |", "|---|---|",
        f"| TM-score | {summary['tmscore_residue_matched']:.3f} |",
        f"| RMSD (optimal superposition, {pair['length_tbm']} residues) | {summary['rmsd_residue_matched_tmscore']:.2f} Å |",
        f"| GDT-TS | {summary['gdt_ts_residue_matched']:.3f} |",
        "",
        "## Structural agreement, sequence-independent (TM-align)",
        "",
        "| Metric | Value |", "|---|---|",
        f"| TM-score (normalized by TBM, {pair['length_tbm']} res) | {pair['tm_score_norm_tbm']:.3f} |",
        f"| TM-score (normalized by AF3, {pair['length_af3']} res) | {pair['tm_score_norm_af3']:.3f} |",
        f"| RMSD over aligned residues | {pair['rmsd']:.2f} Å |",
        f"| Aligned length | {pair['aligned_length']} |",
        f"| Aligned pairs within 5 Å | {pair['pairs_within_5A']} |",
        f"| Aligned pairs at the same residue number | {in_register} of {pair['aligned_length']} |",
        f"| Residue-matched CA RMSD, all {len(dev_by)} residues (same superposition) | {summary['residue_matched_rmsd']:.2f} Å |",
        "",
        f"## Model quality, residues {args.start}-{args.end}",
        "",
        "| Metric | TBM | AlphaFold3 |", "|---|---|---|",
        f"| DOPE (lower = better) | {tbm_dope:.1f} | {af3_dope:.1f} |",
        f"| z-DOPE (lower = better; < -1 native-like) | {tbm_z:.3f} | {af3_z:.3f} |",
        f"| TM-score vs template 4EHX_A | {summary['tm_tbm_vs_template']:.3f} | {summary['tm_af3_vs_template']:.3f} |",
        f"| Mean pLDDT (all / {args.start}-{args.end}) | n/a | {pair['af3_mean_plddt_all']:.1f} / {pair['af3_mean_plddt_aligned']:.1f} |",
        f"| pTM / ranking score | n/a | {top_conf.get('ptm')} / {top_conf.get('ranking_score')} |",
        "",
        "## Regions",
        "",
        "| Region | AF3 pLDDT | TBM DOPE | AF3 DOPE | CA deviation (Å) |", "|---|---|---|---|---|",
    ]
    for r in region_rows:
        md.append(f"| {r['region']} | {r['af3_plddt']} | {r['tbm_dope'] if r['region'] != region_rows[0]['region'] else 'not modeled'} | "
                  f"{r['af3_dope'] if r['region'] != region_rows[0]['region'] else 'not scored'} | "
                  f"{r['ca_deviation'] if r['resolved_pairs'] else 'n/a'} |")
    with open(os.path.join(pres, "summary.md"), "w") as handle:
        handle.write("\n".join(md) + "\n")
    print("\n".join(md))
    print(f"\nWrote figures and tables to {pres}")


if __name__ == "__main__":
    main()
