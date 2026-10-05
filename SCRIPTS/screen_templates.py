"""Screen targets for template quality from their HHsearch results.

Reads targets/<T>/hhsearch_output/<T>.hhr for each target, prints the top
hits, and marks each target PASS or FAIL. A hit qualifies when

  probability >= 95 AND identity >= 30 AND query coverage >= 70

and a target passes when any of its hits qualifies (not only the rank-1 hit).
When several qualify, the best template is the one with the highest identity,
then coverage, then probability. The CSV's probability/E-value/identity/coverage
columns describe the best template for a PASS and the rank-1 hit for a FAIL;
top_hit_pdb is always the rank-1 hit.

Combined coverage is the share of query residues covered by the union of the
query ranges of all hits with probability >= 95. Results are written to a CSV.
"""

import argparse
import csv
import os

from top_hits import parse_hhr

TBM = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def combined_coverage(hits, query_len, min_prob):
    """Percent of query residues covered by any hit with prob >= min_prob, and the hit count."""
    covered = set()
    strong = [h for h in hits if h["prob"] >= min_prob]
    for h in strong:
        covered.update(range(h["q_start"], h["q_end"] + 1))
    return 100 * len(covered) / query_len, len(strong)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("targets", nargs="+")
    parser.add_argument("--targets-dir", default=os.path.join(TBM, "targets"))
    parser.add_argument("--top", type=int, default=3)
    parser.add_argument("--min-prob", type=float, default=95)
    parser.add_argument("--min-identity", type=float, default=30)
    parser.add_argument("--min-coverage", type=float, default=70)
    parser.add_argument("--csv", default=os.path.join(TBM, "targets/screening_summary.csv"))
    args = parser.parse_args()

    rows = []
    for target in args.targets:
        hhr = os.path.join(args.targets_dir, target, "hhsearch_output", f"{target}.hhr")
        query_len, hits = parse_hhr(hhr)

        print(f"\n{target} ({query_len} residues)")
        print(f"{'#':>2}  {'PDB':<8} {'Prob':>6} {'E-value':>9} {'Ident':>6} {'Query':>9} {'Cover':>6}  Protein")
        for h in hits[:args.top]:
            print(f"{h['rank']:>2}  {h['id']:<8} {h['prob']:>6.1f} {h['evalue']:>9.2g} {h['identity']:>5}% "
                  f"{h['q_start']:>4}-{h['q_end']:<4} {h['q_coverage']:>5.0f}%  {h['description'].split(';')[0][:55]}")

        def failures(h):
            return [name for name, ok in (
                ("probability", h["prob"] >= args.min_prob),
                ("identity", h["identity"] >= args.min_identity),
                ("coverage", h["q_coverage"] >= args.min_coverage)) if not ok]

        top = hits[0]
        qualifying = [h for h in hits if not failures(h)]
        best_template = max(qualifying, key=lambda h: (h["identity"], h["q_coverage"], h["prob"]),
                            default=None)
        best = best_template or top
        failed = failures(best)
        combined, n_strong = combined_coverage(hits, query_len, args.min_prob)
        status = "PASS" if best_template else "FAIL"
        if best_template:
            note = f"best template {best['id']} (rank {best['rank']}); {len(qualifying)} qualifying hit(s)"
        else:
            note = f"no qualifying hit; rank-1 hit fails on {', '.join(failed)}"
        if failed == ["coverage"]:
            note += f"; combined coverage of {n_strong} hits with prob >= {args.min_prob:g}: {combined:.0f}%"
        print(f"   -> {status} ({note})")

        rows.append({
            "target": target, "length": query_len, "top_hit_pdb": top["id"],
            "best_template": best_template["id"] if best_template else "",
            "probability": best["prob"], "evalue": best["evalue"], "identity": best["identity"],
            "coverage": round(best["q_coverage"], 1),
            "combined_coverage": round(combined, 1) if n_strong else "",
            "status": status,
        })

    with open(args.csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {args.csv}")


if __name__ == "__main__":
    main()
