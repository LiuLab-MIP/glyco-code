"""Command line: `python -m glyco {extract,analyse,all}`."""

from __future__ import annotations

import argparse
import sys

from . import pipeline
from .cohorts import COHORTS


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m glyco", description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    extract = commands.add_parser(
        "extract", help="read raw images and write the patch tables (needs raw_images/)")
    extract.add_argument("--cohort", action="append", choices=sorted(COHORTS),
                         help="cohort to extract; repeat for several (default: all)")

    commands.add_parser(
        "analyse", help="cross-validate, fit, and write results/ and models/ from the tables in data/")

    everything = commands.add_parser("all", help="extract, then analyse")
    everything.add_argument("--skip-extract", action="store_true",
                            help="use the tables already in data/ instead of re-reading images")

    args = parser.parse_args(argv)
    if args.command == "extract":
        return pipeline.extract(args.cohort)
    if args.command == "analyse":
        return pipeline.analyse()
    return pipeline.run_all(skip_extract=args.skip_extract)


if __name__ == "__main__":
    sys.exit(main())
