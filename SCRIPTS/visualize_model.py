"""Superimpose a MODELLER model on its template and render PNGs with PyMOL.

Writes three files to the output directory:
  <model>_vs_<template>.png  model and template superimposed, one color each
  <model>_dope.png           model colored by smoothed per-residue DOPE
  <model>_dope.profile       MODELLER per-residue DOPE energy profile
and a <model>_dope.pdb copy with DOPE stored in the B-factor column.
"""

import argparse
import os

from modeller import Environ, Selection, log
from modeller.scripts import complete_pdb
from pymol import cmd

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
EXAMPLE = os.path.join(os.path.dirname(SCRIPTS), "examples", "test_target")


def dope_profile(model_pdb, profile_out, window=15):
    """Return per-residue DOPE (normalized, smoothed) in model residue order."""
    log.none()
    env = Environ()
    env.libs.topology.read(file="$(LIB)/top_heav.lib")
    env.libs.parameters.read(file="$(LIB)/par.lib")
    mdl = complete_pdb(env, model_pdb)
    Selection(mdl).assess_dope(output="ENERGY_PROFILE NO_REPORT", file=profile_out,
                               normalize_profile=True, smoothing_window=window)
    values = []
    with open(profile_out) as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                values.append(float(line.split()[-1]))
    return values


def parse_ranges(text):
    """Parse '182-215,331-345' into [(182, 215), (331, 345)]."""
    ranges = []
    for part in text.split(","):
        lo, hi = part.split("-")
        ranges.append((int(lo), int(hi)))
    return ranges


def render(path, width, height):
    cmd.set("ray_opaque_background", 1)
    cmd.png(path, width=width, height=height, dpi=300, ray=1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", default=os.path.join(EXAMPLE, "models/first_pass/test_target.B99990003.pdb"))
    parser.add_argument("--template", default=os.path.join(EXAMPLE, "templates/4ehx_A.pdb"))
    parser.add_argument("--regions", default="182-215,331-345",
                        help="Query residue ranges to outline, e.g. 182-215,331-345")
    parser.add_argument("--window", type=int, default=15, help="DOPE smoothing window")
    parser.add_argument("--out-dir", default=os.path.join(EXAMPLE, "models/first_pass/figures"))
    parser.add_argument("--width", type=int, default=1600)
    parser.add_argument("--height", type=int, default=1200)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(args.model))[0]
    tname = os.path.splitext(os.path.basename(args.template))[0]
    regions = parse_ranges(args.regions)

    # Per-residue DOPE
    profile_file = os.path.join(args.out_dir, f"{stem}_dope.profile")
    dope = dope_profile(args.model, profile_file, args.window)

    # Load and superimpose (sequence-independent structural fit)
    cmd.reinitialize()
    cmd.load(args.model, "target")
    cmd.load(args.template, "template")
    cmd.remove("hydro")
    rms = cmd.super("target and name CA", "template and name CA")

    ca_resv = []
    cmd.iterate("target and name CA", "ca_resv.append(int(resv))", space={"ca_resv": ca_resv})
    if len(ca_resv) != len(dope):
        raise ValueError(f"DOPE profile has {len(dope)} residues but model has {len(ca_resv)}")
    dope_by_resi = dict(zip(ca_resv, dope))
    cmd.alter("target", "b = dope_by_resi[int(resv)]", space={"dope_by_resi": dope_by_resi})
    cmd.save(os.path.join(args.out_dir, f"{stem}_dope.pdb"), "target")

    # Common style
    cmd.bg_color("white")
    cmd.hide("everything")
    cmd.show("cartoon")
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("ray_shadows", 0)
    cmd.set("antialias", 2)
    cmd.orient("target or template")
    cmd.zoom("target or template", 2)

    # Figure 1: model vs template, one color per structure
    cmd.color("marine", "target")
    cmd.color("orange", "template")
    cmd.set("cartoon_transparency", 0.0)
    overlay_png = os.path.join(args.out_dir, f"{stem}_vs_{tname}.png")
    render(overlay_png, args.width, args.height)

    # Figure 2: model colored by DOPE, template as faint reference
    lo, hi = min(dope), max(dope)
    cmd.spectrum("b", "blue_white_red", "target", minimum=lo, maximum=hi)
    cmd.color("grey80", "template")
    cmd.set("cartoon_transparency", 0.75, "template")
    for start, end in regions:
        cmd.label(f"target and resi {(start + end) // 2} and name CA", f"'{start}-{end}'")
    cmd.set("label_size", 28)
    cmd.set("label_color", "black")
    cmd.set("label_outline_color", "white")
    cmd.set("float_labels", 1)
    dope_png = os.path.join(args.out_dir, f"{stem}_dope.png")
    render(dope_png, args.width, args.height)

    # Report
    mean = sum(dope) / len(dope)
    print(f"Superposition: RMSD {rms[0]:.2f} A over {rms[1]} CA atoms "
          f"(after refinement; {rms[4]} before outlier rejection), template {tname}")
    print(f"Per-residue DOPE (window {args.window}): min {lo:.4f}  max {hi:.4f}  mean {mean:.4f}")
    for start, end in regions:
        vals = [dope_by_resi[r] for r in range(start, end + 1) if r in dope_by_resi]
        print(f"  {start}-{end}: mean {sum(vals) / len(vals):.4f}  max {max(vals):.4f}")
    worst = sorted(dope_by_resi.items(), key=lambda kv: kv[1], reverse=True)[:10]
    print("  Highest-DOPE residues:", ", ".join(f"{r} ({v:.3f})" for r, v in worst))
    print(f"Wrote {overlay_png}")
    print(f"Wrote {dope_png}")


if __name__ == "__main__":
    main()
