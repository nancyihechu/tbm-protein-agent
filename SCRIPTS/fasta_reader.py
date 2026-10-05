"""Read protein sequences from FASTA files."""

import sys


def read_fasta(path):
    """Return the sequence from a single-record FASTA file.

    Header lines (starting with '>') and blank lines are skipped, and
    sequence lines are joined into one uppercase string with whitespace
    removed. Raises ValueError if the file has no sequence or contains
    more than one record.
    """
    records = read_fasta_records(path)
    if not records:
        raise ValueError(f"No sequence found in {path}")
    if len(records) > 1:
        raise ValueError(
            f"{path} contains {len(records)} records; use read_fasta_records()"
        )
    return records[0][1]


def read_fasta_records(path):
    """Return a list of (header, sequence) tuples from a FASTA file."""
    records = []
    header = None
    chunks = []

    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(";"):
                continue
            if line.startswith(">"):
                if header is not None or chunks:
                    records.append((header, "".join(chunks).upper()))
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append("".join(line.split()))

    if header is not None or chunks:
        records.append((header, "".join(chunks).upper()))

    # Drop records with a header but no sequence
    return [(h, s) for h, s in records if s]


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: python fasta_reader.py <file.fasta>")
    for header, seq in read_fasta_records(sys.argv[1]):
        print(f">{header}")
        print(f"Length: {len(seq)}")
        print(seq)
