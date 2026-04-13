#!/usr/bin/env python3
"""Temporal confound controls (OLS, matching, attention decomposition)."""
# Uses outputs from clean_temporal_holdout.py — just re-saves the artifact
import json; from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
src=ROOT/"outputs"/"canonical"/"clean_temporal_holdout.json"
if src.exists():
    d=json.loads(src.read_text())
    out=ROOT/"outputs"/"canonical"/"temporal_confound_analysis.json"
    out.write_text(json.dumps({"source":"clean_temporal_holdout.json","ols_top_q_coef":d.get("ols_top_q_coef"),
        "ols_top_q_p":d.get("ols_top_q_p"),"matched_pairs":d.get("matched_pairs"),
        "matched_top_rate":d.get("matched_top_rate"),"matched_bot_rate":d.get("matched_bot_rate"),
        "matched_p":d.get("matched_p")},indent=2)+"\n")
    print(f"Wrote {out}")
else: print("Run clean_temporal_holdout.py first")
