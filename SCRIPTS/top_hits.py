"""Report the top hits of an HHsearch .hhr file.

For each hit: template ID, probability, E-value, sequence identity, aligned
columns, and the query and template residue ranges, with query coverage as a
percentage of the query length.
"""

import argparse
import csv
import re
import sys

ALN_LINE = re.compile(r"^([QT]) (\S+)\s+(\d+) (\S+)\s+(\d+) \((\d+)\)")


def parse_hhr(path):
    """Return (query_length, hits) where hits is a list of dicts in rank order."""
    with open(path) as handle:
        lines = handle.read().splitlines()

    query_len = int(next(l.split()[1] for l in lines if l.startswith("Match_columns")))
    hits, hit = [], None
    for line in lines:
        if re.match(r"^No \d+$", line):
            hit = {"rank": int(line.split()[1]), "q": [], "t": []}
            hits.append(hit)
        elif hit is None:
            continue
        elif line.startswith(">") and "id" not in hit:
            parts = line[1:].split(None, 1)
            hit["id"] = parts[0]
            hit["description"] = parts[1] if len(parts) > 1 else ""
        elif line.startswith("Probab="):
            fields = dict(f.split("=") for f in line.split())
            hit["prob"] = float(fields["Probab"])
            hit["evalue"] = float(fields["E-value"])
            hit["identity"] = int(fields["Identities"].rstrip("%"))
            hit["aligned_cols"] = int(fields["Aligned_cols"])
        else:
            m = ALN_LINE.match(line)
            if m and m.group(2) not in ("Consensus", "ss_dssp", "ss_pred"):
                rng = (int(m.group(3)), int(m.group(5)), int(m.group(6)))
                hit["q" if m.group(1) == "Q" else "t"].append(rng)

    for hit in hits:
        hit["q_start"], hit["q_end"] = hit["q"][0][0], hit["q"][-1][1]
        hit["t_start"], hit["t_end"], hit["t_len"] = hit["t"][0][0], hit["t"][-1][1], hit["t"][0][2]
        hit["q_coverage"] = 100 * (hit["q_end"] - hit["q_start"] + 1) / query_len
    return query_len, hits


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("hhr", nargs="+")
    parser.add_argument("-n", "--top", type=int, default=3)
    parser.add_argument("--csv", help="Also write all reported hits to this CSV file")
    args = parser.parse_args()

    rows = []
    for path in args.hhr:
        query_len, hits = parse_hhr(path)
        query = next(l.split()[1] for l in open(path) if l.startswith("Query"))
        print(f"\n{query} ({query_len} residues)  {path}")
        print(f"{'#':>2}  {'Template':<8} {'Prob':>6} {'E-value':>9} {'Ident':>5} {'Cols':>5}  "
              f"{'Query':>9} {'Cover':>6}  {'Template':>12}  Description")
        for h in hits[:args.top]:
            print(f"{h['rank']:>2}  {h['id']:<8} {h['prob']:>6.1f} {h['evalue']:>9.2g} {h['identity']:>4}% "
                  f"{h['aligned_cols']:>5}  {h['q_start']:>4}-{h['q_end']:<4} {h['q_coverage']:>5.0f}%  "
                  f"{h['t_start']:>4}-{h['t_end']:<4}({h['t_len']})  {h['description'][:60]}")
            rows.append({"query": query, "query_len": query_len, **{k: h[k] for k in (
                "rank", "id", "prob", "evalue", "identity", "aligned_cols", "q_start", "q_end",
                "t_start", "t_end", "t_len", "description")}, "q_coverage": round(h["q_coverage"], 1)})

    if args.csv:
        with open(args.csv, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


if __name__ == "__main__":
    sys.exit(main())
