#!/usr/bin/env python3
"""Within-compound dose-response analysis for top-ranked robust compounds."""
from __future__ import annotations

import csv
import json
import re
from collections import defaultdict
from pathlib import Path

from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]


def parse_dose(dose_str: str) -> tuple[float, str] | None:
    dose_str = dose_str.strip()
    m = re.match(
        r"^([\d.]+)\s*"
        r"(mM|µM|μM|uM|nM|ppm|%|mg/kg|mg/ml|mg/mL|mg/l|mg/L|"
        r"µg/ml|µg/mL|μg/ml|μg/mL|ug/ml|mg|µg|ug|g|g/l)",
        dose_str,
        re.I,
    )
    if not m:
        return None
    val = float(m.group(1))
    unit = m.group(2).lower().replace("μ", "µ")
    return (val, unit)


def main() -> int:
    drugage_path = ROOT / "data" / "drugage_build5_2024-11-29.csv"
    rankings_path = ROOT / "outputs" / "canonical" / "robustness_rankings.csv"
    output_path = ROOT / "outputs" / "canonical" / "dose_response_analysis.json"

    rows: list[dict[str, str]] = []
    with open(drugage_path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rows.append(row)

    rankings: list[dict[str, str]] = []
    with open(rankings_path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rankings.append(row)

    top_n = 20
    top_compounds = [
        r.get("compound", r.get("compound_name", "")).strip()
        for r in rankings[:top_n]
    ]

    results: dict[str, dict] = {}
    for compound in top_compounds:
        compound_rows = [r for r in rows if r["compound_name"].strip() == compound]
        pairs: list[dict] = []
        for r in compound_rows:
            parsed = parse_dose(r.get("dosage", ""))
            if parsed is None:
                continue
            try:
                effect = float(r["avg_lifespan_change_percent"])
            except (ValueError, KeyError):
                continue
            pairs.append({
                "dose_val": parsed[0],
                "dose_unit": parsed[1],
                "effect": effect,
                "species": r["species"],
            })

        if len(pairs) < 3:
            results[compound] = {
                "n_parseable": len(pairs),
                "n_total": len(compound_rows),
                "status": "insufficient_data",
            }
            continue

        by_unit: dict[str, list] = defaultdict(list)
        for p in pairs:
            by_unit[p["dose_unit"]].append(p)

        dominant_unit = max(by_unit.keys(), key=lambda u: len(by_unit[u]))
        unit_pairs = by_unit[dominant_unit]

        if len(unit_pairs) < 3:
            results[compound] = {
                "n_parseable": len(pairs),
                "n_total": len(compound_rows),
                "dominant_unit": dominant_unit,
                "n_same_unit": len(unit_pairs),
                "status": "insufficient_same_unit",
            }
            continue

        doses = [p["dose_val"] for p in unit_pairs]
        effects = [p["effect"] for p in unit_pairs]

        if len(set(doses)) < 2:
            results[compound] = {
                "n_parseable": len(pairs),
                "n_total": len(compound_rows),
                "dominant_unit": dominant_unit,
                "n_same_unit": len(unit_pairs),
                "status": "constant_dose",
            }
            continue

        rho, pval = spearmanr(doses, effects)
        results[compound] = {
            "n_parseable": len(pairs),
            "n_total": len(compound_rows),
            "dominant_unit": dominant_unit,
            "n_same_unit": len(unit_pairs),
            "n_units": len(by_unit),
            "units": sorted(by_unit.keys()),
            "spearman_rho": round(float(rho), 4) if rho == rho else None,
            "spearman_p": round(float(pval), 4) if pval == pval else None,
            "dose_min": min(doses),
            "dose_max": max(doses),
            "effect_min": round(min(effects), 2),
            "effect_max": round(max(effects), 2),
            "status": "computed",
        }

    computed = [r for r in results.values() if r["status"] == "computed"]
    valid_rho = [r for r in computed if r["spearman_rho"] is not None]

    summary = {
        "top_n": top_n,
        "computed": len(computed),
        "constant_dose": sum(1 for r in results.values() if r["status"] == "constant_dose"),
        "insufficient": sum(1 for r in results.values() if r["status"] in ("insufficient_data", "insufficient_same_unit")),
        "positive_correlation": sum(1 for r in valid_rho if r["spearman_rho"] > 0),
        "negative_correlation": sum(1 for r in valid_rho if r["spearman_rho"] < 0),
        "sig_positive_005": sum(1 for r in valid_rho if r["spearman_rho"] > 0 and r["spearman_p"] < 0.05),
        "sig_negative_005": sum(1 for r in valid_rho if r["spearman_rho"] < 0 and r["spearman_p"] < 0.05),
    }

    output = {"summary": summary, "compounds": results}
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {output_path}")
    print(f"Computed: {summary['computed']}/{top_n}")
    print(f"Positive dose-response: {summary['positive_correlation']}, Negative: {summary['negative_correlation']}")
    print(f"Significant positive (p<0.05): {summary['sig_positive_005']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
