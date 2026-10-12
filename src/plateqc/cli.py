"""plateqc command line interface."""

import argparse
import json
import sys

from .qc import load_plate, run_qc, qc_report, load_layout


def cmd_qc(args):
    try:
        layout = load_layout(args.layout)
    except (OSError, ValueError) as e:
        print(f"plateqc: {e}", file=sys.stderr)
        return 2
    try:
        wells, info = load_plate(args.plate, value_column=args.value_column)
    except (OSError, ValueError) as e:
        print(f"plateqc: {e}", file=sys.stderr)
        return 2
    report = run_qc(wells, layout, min_controls=args.min_controls)
    report["plate"] = args.plate
    report["plate_info"] = {k: v for k, v in info.items() if k != "warnings"}
    if args.format == "json":
        print(json.dumps(report, indent=2, default=str))
    else:
        sys.stdout.write(f"plate: {args.plate} ({info['format']}, {info['n_wells']} wells)\n")
        sys.stdout.write(qc_report(report))
    bad = any("poor" in w for w in report["warnings"])
    assay = report.get("assay", {})
    fail = args.fail_below is not None and (
        assay.get("zprime") is None or assay["zprime"] < args.fail_below
    )
    return 1 if (bad or fail) else 0


def build_parser():
    p = argparse.ArgumentParser(
        prog="plateqc",
        description="Microplate assay QC: Z'-factor, SSMD, and spatial diagnostics.",
    )
    sub = p.add_subparsers(dest="cmd", required=True)
    q = sub.add_parser("qc", help="run QC on a plate reader CSV")
    q.add_argument("plate", help="plate reader export CSV (grid or tidy format)")
    q.add_argument("--layout", required=True,
                   help="JSON file mapping control groups to well/range specs")
    q.add_argument("--value-column", default=None,
                   help="tidy-format value column (auto-detected if omitted)")
    q.add_argument("--min-controls", type=int, default=4,
                   help="warn below this many control values (default 4)")
    q.add_argument("--format", choices=("text", "json"), default="text")
    q.add_argument("--fail-below", type=float, default=None,
                   help="exit 1 if Z'-factor is below this threshold")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.cmd == "qc":
        return cmd_qc(args)
    return 2


if __name__ == "__main__":
    sys.exit(main())
