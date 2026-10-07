"""Run from the repository root: python -m agent2 audit|evaluate."""
import argparse
from pathlib import Path

from .workflow import run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["audit", "evaluate"])
    parser.add_argument("--config", default="agent2/config.json", help="Repository-relative JSON config")
    parser.add_argument("--out", help="New directory below agent2/results; existing directories are refused")
    parser.add_argument("--tmscore", help="Official TMscore executable; not TMalign")
    parser.add_argument("--require-complete", action="store_true", help="Exit 2 if any assignment deliverable is pending")
    args = parser.parse_args()
    try:
        summary, output = run(Path(__file__).resolve().parent.parent, args)
    except (ValueError, OSError, KeyError) as exc:
        parser.exit(1, f"Agent 2: {exc}\n")
    print(f"Agent 2: {summary['status']}\nReport: {output / 'REPORT.md'}")
    if summary.get("errors"):
        return 1
    if args.require_complete and summary["status"] != "COMPLETE":
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

