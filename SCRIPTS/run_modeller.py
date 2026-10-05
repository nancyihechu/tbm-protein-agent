"""Build homology models with MODELLER automodel from a PIR alignment.

The alignment is first trimmed to a query residue range (columns before the
start or after the end of the range are dropped, and the template's
structureX start/end residues are adjusted to match). Models are renumbered
so residue numbers match the full query sequence, then scored with DOPE and
GA341.
"""

import argparse
import os
import re
import sys

from modeller import Environ, log, parallel
from modeller.automodel import AutoModel, assess

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
EXAMPLE = os.path.join(os.path.dirname(SCRIPTS), "examples", "test_target")


def read_pir(path):
    """Return a list of [code, header, sequence] entries from a PIR file."""
    entries = []
    with open(path) as handle:
        blocks = handle.read().split(">P1;")[1:]
    for block in blocks:
        lines = block.strip().splitlines()
        entries.append([lines[0].strip(), lines[1].strip(),
                        "".join(lines[2:]).replace(" ", "").rstrip("*")])
    return entries


def pdb_residue_ids(pdb_file, chain):
    """Return residue numbers (with insertion codes) of a chain, in file order."""
    ids = []
    with open(pdb_file) as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")) and line[21] == chain:
                rid = line[22:27].strip()
                if not ids or ids[-1] != rid:
                    ids.append(rid)
    return ids


def trim_alignment(pir_in, pir_out, query_code, start, end, template_dir):
    """Write a PIR limited to query residues start..end (1-based, inclusive)."""
    entries = read_pir(pir_in)
    query = next(e for e in entries if e[0] == query_code)

    # Columns covering query residues start..end
    pos, keep_from, keep_to = 0, None, None
    for col, char in enumerate(query[2]):
        if char != "-":
            pos += 1
            if pos == start and keep_from is None:
                keep_from = col
            if pos == end:
                keep_to = col
    if keep_from is None or keep_to is None:
        raise ValueError(f"Query has fewer than {end} residues")

    # Keep template-only overhang columns after the last query residue
    if end == pos:
        keep_to = len(query[2]) - 1

    out = []
    for code, header, seq in entries:
        fields = header.split(":")
        if fields[0].startswith("structure"):
            dropped_before = len(seq[:keep_from].replace("-", ""))
            kept = len(seq[keep_from:keep_to + 1].replace("-", ""))
            ids = pdb_residue_ids(os.path.join(template_dir, f"{fields[1]}.pdb"), fields[3])
            offset = ids.index(fields[2]) if fields[2] not in ("", "FIRST") else 0
            fields[2] = ids[offset + dropped_before]
            fields[4] = ids[offset + dropped_before + kept - 1]
        out.append((code, ":".join(fields), seq[keep_from:keep_to + 1]))

    # Drop columns that became all-gap
    cols = [i for i in range(len(out[0][2])) if any(e[2][i] != "-" for e in out)]
    with open(pir_out, "w") as handle:
        for code, header, seq in out:
            trimmed = "".join(seq[i] for i in cols)
            wrapped = "\n".join(trimmed[i:i + 75] for i in range(0, len(trimmed), 75))
            handle.write(f">P1;{code}\n{header}\n{wrapped}*\n\n")


class RenumberedModel(AutoModel):
    first_residue = 1

    def special_patches(self, aln):
        self.rename_segments(segment_ids=["A"], renumber_residues=[self.first_residue])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment", default=os.path.join(EXAMPLE, "alignments/test_target_4ehx_A.pir"))
    parser.add_argument("--template", nargs="+", default=["4ehx_A"],
                        help="Template code(s) as in the PIR file, e.g. 4ehx_A")
    parser.add_argument("--template-dir", default=os.path.join(EXAMPLE, "templates"))
    parser.add_argument("--query", default="test_target")
    parser.add_argument("--start", type=int, default=103)
    parser.add_argument("--end", type=int, default=460)
    parser.add_argument("--n-models", type=int, default=5)
    parser.add_argument("--cpus", type=int, default=1)
    parser.add_argument("--out-dir", default=os.path.join(EXAMPLE, "models/first_pass"))
    args = parser.parse_args()

    template_dir = os.path.abspath(args.template_dir)
    alignment = os.path.abspath(args.alignment)
    os.makedirs(args.out_dir, exist_ok=True)
    os.chdir(args.out_dir)

    trimmed = f"{args.query}_{args.start}-{args.end}_{'+'.join(args.template)}.pir"
    trim_alignment(alignment, trimmed, args.query, args.start, args.end, template_dir)

    log.verbose()
    env = Environ()
    env.io.atom_files_directory = [template_dir]

    model = RenumberedModel(env, alnfile=trimmed, knowns=tuple(args.template),
                            sequence=args.query,
                            assess_methods=[assess.DOPE, assess.GA341])
    model.first_residue = args.start
    model.starting_model = 1
    model.ending_model = args.n_models

    if args.cpus > 1:
        job = parallel.Job()
        for _ in range(args.cpus):
            job.append(parallel.LocalWorker())
        model.use_parallel_job(job)

    model.make()

    ok = [m for m in model.outputs if m["failure"] is None]
    failed = [m for m in model.outputs if m["failure"] is not None]
    print("\n=== Summary ===")
    print(f"{'Model':<30}{'DOPE':>14}{'GA341':>10}")
    for m in ok:
        print(f"{m['name']:<30}{m['DOPE score']:>14.3f}{m['GA341 score'][0]:>10.3f}")
    for m in failed:
        print(f"{m['name']:<30}  FAILED: {m['failure']}")
    if ok:
        best = min(ok, key=lambda m: m["DOPE score"])
        print(f"Best DOPE: {best['name']} ({best['DOPE score']:.3f})")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
