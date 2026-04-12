"""Command-line interface."""

from __future__ import annotations

import argparse

from .anage import run_anage_context_report
from .config import load_config
from .constants import DEFAULT_CONFIG
from .pipeline import run_pipeline
from .utils import set_runtime_environment
from .verify import verify_run


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Claim-certified DrugAge robustness skill")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="Execute the canonical DrugAge pipeline")
    run_parser.add_argument("--config", default=DEFAULT_CONFIG)
    run_parser.add_argument("--out", required=True)
    run_parser.add_argument("--gender", choices=["Male", "Female", "Hermaphrodite", "Pooled", "Unknown"],
                           help="Filter experiments to a single gender before scoring")

    verify_parser = subparsers.add_parser("verify", help="Verify a canonical run directory")
    verify_parser.add_argument("--config", default=DEFAULT_CONFIG)
    verify_parser.add_argument("--run-dir", required=True)

    anage_parser = subparsers.add_parser("anage-context-report", help="Generate the optional descriptive AnAge report")
    anage_parser.add_argument("--config", default=DEFAULT_CONFIG)
    anage_parser.add_argument("--run-dir", required=True)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    config = load_config(args.config)
    set_runtime_environment(int(config.runtime["pythonhashseed"]))
    if args.command == "run":
        if hasattr(args, 'gender') and args.gender:
            config.raw["pipeline"]["gender_filter"] = args.gender
        run_pipeline(config, args.out)
        return 0
    if args.command == "verify":
        verification = verify_run(config, args.run_dir)
        return 0 if verification["status"] == "passed" else 1
    if args.command == "anage-context-report":
        run_anage_context_report(config, args.run_dir)
        return 0
    parser.error("Unknown command")
    return 2
