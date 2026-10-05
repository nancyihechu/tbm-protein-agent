"""Build a MODELLER PIR alignment from one or more hits in an HHsearch .hhr file.

Each query-template pairing is taken from the hit's alignment block. Template
positions are HHsearch SEQRES indices, so they are mapped onto residues that
actually have coordinates in the mmCIF; unobserved template residues become
gaps. With several hits, the pairwise alignments are merged with the query as
the anchor: template residues aligned to the same query residue share a
column, and template insertions get their own columns. Each template keeps
only its aligned range; the query is always complete.
"""

import argparse
import os
import re

from fasta_reader import read_fasta
from prepare_templates import cif_path, read_chain_residues

ALN_LINE = re.compile(r"^([QT]) (\S+)\s+(\d+) (\S+)\s+(\d+) \((\d+)\)")


def parse_hit_alignment(hhr_file, hit_no):
    """Return (template_name, q_start, q_aln, t_start, t_aln) for a hit."""
    with open(hhr_file) as handle:
        lines = handle.read().splitlines()

    try:
        start = lines.index(f"No {hit_no}")
    except ValueError:
        raise ValueError(f"Hit {hit_no} not found in {hhr_file}")

    name = lines[start + 1][1:].split()[0]
    q_parts, t_parts = [], []
    q_start = t_start = None
    for line in lines[start + 2:]:
        if line.startswith("No ") or line.startswith("Done"):
            break
        m = ALN_LINE.match(line)
        if not m:
            continue
        kind, seq_name, begin, aln = m.group(1), m.group(2), int(m.group(3)), m.group(4)
        if seq_name in ("Consensus", "ss_dssp", "ss_pred"):
            continue
        if kind == "Q":
            q_start = begin if q_start is None else q_start
            q_parts.append(aln)
        else:
            t_start = begin if t_start is None else t_start
            t_parts.append(aln)

    q_aln, t_aln = "".join(q_parts), "".join(t_parts)
    if not q_aln or len(q_aln) != len(t_aln):
        raise ValueError(f"Could not parse alignment for hit {hit_no}")
    return name, q_start, q_aln, t_start, t_aln


def build_columns(query, q_start, q_aln, seqres, t_start, t_aln):
    """Return alignment columns as (query_pos, template_pos) pairs (None = gap)."""
    t_ungapped = t_aln.replace("-", "")
    expected = seqres[t_start - 1:t_start - 1 + len(t_ungapped)]
    if t_ungapped != expected:
        raise ValueError("HHsearch template sequence does not match mmCIF SEQRES")
    q_ungapped = q_aln.replace("-", "")
    if q_ungapped != query[q_start - 1:q_start - 1 + len(q_ungapped)]:
        raise ValueError("HHsearch query sequence does not match the FASTA")

    cols = [(None, t) for t in range(1, t_start)]
    cols += [(q, None) for q in range(1, q_start)]
    qi, ti = q_start, t_start
    for qc, tc in zip(q_aln, t_aln):
        cols.append((qi if qc != "-" else None, ti if tc != "-" else None))
        qi += qc != "-"
        ti += tc != "-"
    cols += [(q, None) for q in range(qi, len(query) + 1)]
    cols += [(None, t) for t in range(ti, len(seqres) + 1)]
    return cols


def wrap(seq, width=75):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def template_track(query, q_start, q_aln, seqres, t_start, t_aln, observed):
    """Map one pairwise alignment onto the query.

    Returns (aligned, inserts, t_range): aligned[q] is the template residue at
    query position q, inserts[q] the template residues inserted after query
    position q, and t_range the SEQRES positions of observed aligned-range residues.
    """
    build_columns(query, q_start, q_aln, seqres, t_start, t_aln)  # sequence checks
    aligned, inserts, used = {}, {}, []
    qi, ti = q_start, t_start
    for qc, tc in zip(q_aln, t_aln):
        if tc != "-" and ti in observed:
            used.append(ti)
            if qc != "-":
                aligned[qi] = observed[ti][2]
            else:
                inserts.setdefault(qi - 1, []).append(observed[ti][2])
        qi += qc != "-"
        ti += tc != "-"
    return aligned, inserts, used


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("hhr")
    parser.add_argument("query_fasta")
    parser.add_argument("--hit", type=int, nargs="+", default=[1],
                        help="Hit number(s) from the .hhr file (the 'No N' blocks)")
    parser.add_argument("--query-name", default=None)
    parser.add_argument("-o", "--out", default=None)
    args = parser.parse_args()

    query = read_fasta(args.query_fasta)
    qname = args.query_name or os.path.splitext(os.path.basename(args.query_fasta))[0]

    tracks = []
    for hit in args.hit:
        name, q_start, q_aln, t_start, t_aln = parse_hit_alignment(args.hhr, hit)
        pdb_id, chain = name.split("_")
        seqres, observed = read_chain_residues(cif_path(pdb_id), chain)
        aligned, inserts, used = template_track(query, q_start, q_aln, seqres, t_start, t_aln, observed)
        if not used:
            raise ValueError(f"Hit {hit} ({name}) has no aligned residues with coordinates")
        first, last = used[0], used[-1]
        tracks.append({"hit": hit, "name": name, "code": f"{pdb_id.lower()}_{chain}", "chain": chain,
                       "aligned": aligned, "inserts": inserts,
                       "first": f"{observed[first][0]}{observed[first][1]}",
                       "last": f"{observed[last][0]}{observed[last][1]}",
                       "expected": "".join(observed[t][2] for t in sorted(observed) if first <= t <= last),
                       "q_range": (q_start, q_start + len(q_aln.replace("-", "")) - 1),
                       "t_range": (t_start, t_start + len(t_aln.replace("-", "")) - 1)})

    # Query-anchored merge: insertion block after each query position, then the next residue
    rows = {tr["code"]: [] for tr in tracks}
    q_row = []
    for q in range(0, len(query) + 1):
        if q > 0:
            q_row.append(query[q - 1])
            for tr in tracks:
                rows[tr["code"]].append(tr["aligned"].get(q, "-"))
        width = max(len(tr["inserts"].get(q, [])) for tr in tracks)
        if width:
            q_row.extend("-" * width)
            for tr in tracks:
                ins = tr["inserts"].get(q, [])
                rows[tr["code"]].extend(ins + ["-"] * (width - len(ins)))

    blocks = []
    for tr in tracks:
        row = "".join(rows[tr["code"]])
        if row.replace("-", "") != tr["expected"]:
            raise ValueError(f"{tr['code']}: template row does not match its observed residues")
        blocks.append(f">P1;{tr['code']}\n"
                      f"structureX:{tr['code']}:{tr['first']}:{tr['chain']}:{tr['last']}:{tr['chain']}:::-1.00:-1.00\n"
                      f"{wrap(row)}*\n")
    blocks.append(f">P1;{qname}\nsequence:{qname}:::::::0.00: 0.00\n{wrap(''.join(q_row))}*\n")

    out = args.out or f"{qname}_{'+'.join(tr['code'] for tr in tracks)}.pir"
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, "w") as handle:
        handle.write("\n".join(blocks))

    for tr in tracks:
        row = rows[tr["code"]]
        pairs = [(a, b) for a, b in zip(q_row, row) if a != "-" and b != "-"]
        ident = sum(a == b for a, b in pairs)
        print(f"Hit {tr['hit']}: {tr['name']}  query {tr['q_range'][0]}-{tr['q_range'][1]}  "
              f"template SEQRES {tr['t_range'][0]}-{tr['t_range'][1]}  residues {tr['first']}-{tr['last']}  "
              f"aligned pairs {len(pairs)}  identical {ident} ({100 * ident / len(pairs):.1f}%)")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
