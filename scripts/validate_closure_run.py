#!/usr/bin/env python3
"""Validate one HS_CTC closure-event shard and create its success marker."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from closure_common import load_run, validate_run, write_json


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--elastic-tolerance", type=float, default=5.0e-3)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--mark-success", action="store_true")
    args = parser.parse_args()

    try:
        run = load_run(args.run_dir)
        report = validate_run(run, elastic_tolerance=args.elastic_tolerance)
    except Exception as exc:  # validation CLI must produce a machine-readable failure
        report = {"status": "fail", "errors": [str(exc)], "warnings": [], "n_events": 0}

    report_path = args.report or args.run_dir / "validation.json"
    write_json(report_path, report)
    if report["status"] == "pass" and args.mark_success:
        (args.run_dir / "_SUCCESS").write_text("validated\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
