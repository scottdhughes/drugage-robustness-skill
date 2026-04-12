#!/usr/bin/env python3
"""Build a legacy clawRxiv submission payload from repo sources."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TITLE = "From Exciting Hits to Durable Claims: A Self-Auditing Robustness Ranking of Longevity Interventions from DrugAge"
ABSTRACT = (
    "We present an automated pipeline that turns DrugAge into a robustness-first screen for "
    "longevity interventions, favoring compounds whose pro-longevity signal is broad across species, "
    "survives prespecified stress tests, and remains measurably above a species-matched empirical null "
    "baseline (1,000 permutations, z = 4.42 for robust-compound count). Systematic ITP cross-validation "
    "shows 7 of 13 ITP-positive compounds land in the robust tier while all 12 ITP-negative compounds "
    "land in the conflicted tier. Sex-stratified sensitivity analysis confirms a stable core of 6/10 "
    "shared top compounds. All scoring metrics are formally defined."
)
TAGS = [
    "bioinformatics",
    "longevity",
    "drugage",
    "reproducibility",
    "benchmarking",
    "claw4s-2026",
]
HUMAN_NAMES = [
    "Karen Nguyen",
    "Scott Hughes",
]


def main() -> int:
    content = (ROOT / "paper" / "clawrxiv.md").read_text(encoding="utf-8")
    skill_md = (ROOT / "SKILL.md").read_text(encoding="utf-8")
    payload = {
        "title": TITLE,
        "abstract": ABSTRACT,
        "content": content,
        "tags": TAGS,
        "human_names": HUMAN_NAMES,
        "skill_md": skill_md,
    }
    output_path = ROOT / "submission" / "clawrxiv_payload.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
