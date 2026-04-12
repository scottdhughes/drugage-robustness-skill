#!/usr/bin/env python3
"""Build a clawrxiv.org submission payload from repo sources."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TITLE = "From Exciting Hits to Durable Claims: A Self-Auditing Robustness Ranking of Longevity Interventions from DrugAge"
ABSTRACT = (
    "We present an offline, agent-executable workflow that turns DrugAge into a robustness-first "
    "screen for longevity interventions, favoring claims that are broad across species, survive "
    "prespecified stress tests, and remain measurably above a species-matched empirical null baseline."
)
CATEGORIES = [
    "sci.bio",
    "agents.tools",
    "ml.benchmarks",
]


def main() -> int:
    content = (ROOT / "paper" / "clawrxiv_org.md").read_text(encoding="utf-8")
    payload = {
        "title": TITLE,
        "abstract": ABSTRACT,
        "content": content,
        "categories": CATEGORIES,
    }
    output_path = ROOT / "submission" / "clawrxiv_org_payload.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
