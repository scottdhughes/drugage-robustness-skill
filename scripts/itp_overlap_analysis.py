#!/usr/bin/env python3
"""Systematic ITP overlap analysis against DrugAge robustness tiers."""
from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    drugage_path = ROOT / "data" / "drugage_build5_2024-11-29.csv"
    rankings_path = ROOT / "outputs" / "canonical" / "robustness_rankings.csv"
    output_path = ROOT / "outputs" / "canonical" / "itp_overlap_analysis.json"

    rows: list[dict[str, str]] = []
    with open(drugage_path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rows.append(row)

    rankings: list[dict[str, str]] = []
    with open(rankings_path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rankings.append(row)

    compound_tier: dict[str, str] = {}
    compound_rank: dict[str, int] = {}
    compound_score: dict[str, float] = {}
    for i, r in enumerate(rankings):
        name = r.get("compound", r.get("compound_name", "")).strip()
        compound_tier[name] = r.get("tier", r.get("evidence_tier", "")).strip()
        compound_rank[name] = i + 1
        compound_score[name] = float(r.get("robustness_score", 0))

    itp_compounds = set(
        r["compound_name"].strip()
        for r in rows
        if r.get("ITP", "").strip() == "Yes"
    )

    itp_positive: set[str] = set()
    itp_negative: set[str] = set()
    itp_mixed: set[str] = set()
    for c in itp_compounds:
        itp_rows = [
            r for r in rows
            if r["compound_name"].strip() == c and r.get("ITP", "").strip() == "Yes"
        ]
        effects: list[float] = []
        for r in itp_rows:
            try:
                effects.append(float(r["avg_lifespan_change_percent"]))
            except (ValueError, KeyError):
                pass
        if not effects:
            continue
        pos = sum(1 for e in effects if e > 0)
        neg = sum(1 for e in effects if e <= 0)
        if pos > 0 and neg == 0:
            itp_positive.add(c)
        elif neg > 0 and pos == 0:
            itp_negative.add(c)
        else:
            itp_mixed.add(c)

    top_quartile_cutoff = len(rankings) // 4
    itp_matched = [c for c in itp_compounds if c in compound_tier]

    by_tier = Counter(compound_tier[c] for c in itp_matched)
    in_top_quartile = [c for c in itp_matched if compound_rank.get(c, 9999) <= top_quartile_cutoff]

    positive_by_tier = Counter(compound_tier.get(c, "?") for c in itp_positive if c in compound_tier)
    negative_by_tier = Counter(compound_tier.get(c, "?") for c in itp_negative if c in compound_tier)

    compound_details = sorted(
        [
            {
                "compound": c,
                "rank": compound_rank.get(c),
                "tier": compound_tier.get(c),
                "robustness_score": compound_score.get(c),
                "itp_direction": (
                    "positive" if c in itp_positive
                    else "negative" if c in itp_negative
                    else "mixed"
                ),
            }
            for c in itp_matched
        ],
        key=lambda x: x["rank"] or 9999,
    )

    result = {
        "total_itp_compounds_in_drugage": len(itp_compounds),
        "itp_matched_in_rankings": len(itp_matched),
        "itp_by_tier": dict(by_tier),
        "itp_positive_count": len(itp_positive),
        "itp_negative_count": len(itp_negative),
        "itp_mixed_count": len(itp_mixed),
        "itp_positive_by_tier": dict(positive_by_tier),
        "itp_negative_by_tier": dict(negative_by_tier),
        "top_quartile_cutoff_rank": top_quartile_cutoff,
        "itp_in_top_quartile": len(in_top_quartile),
        "itp_in_top_quartile_fraction": round(len(in_top_quartile) / max(len(itp_matched), 1), 4),
        "itp_positive_in_robust_count": positive_by_tier.get("robust", 0),
        "itp_positive_in_robust_fraction": round(
            positive_by_tier.get("robust", 0) / max(len(itp_positive), 1), 4
        ),
        "itp_negative_all_conflicted": negative_by_tier.get("conflicted", 0) == len(itp_negative),
        "compound_details": compound_details,
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {output_path}")
    print(f"ITP compounds matched: {len(itp_matched)}/{len(itp_compounds)}")
    print(f"ITP-positive in robust tier: {positive_by_tier.get('robust', 0)}/{len(itp_positive)} ({result['itp_positive_in_robust_fraction']:.1%})")
    print(f"ITP-negative all conflicted: {result['itp_negative_all_conflicted']}")
    print(f"ITP in top quartile: {len(in_top_quartile)}/{len(itp_matched)} ({result['itp_in_top_quartile_fraction']:.1%})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
