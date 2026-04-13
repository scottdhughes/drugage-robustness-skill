#!/usr/bin/env python3
"""Compound class (pure vs extract) sensitivity."""
import csv,json,re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
rankings=list(csv.DictReader(open(ROOT/"outputs"/"canonical"/"robustness_rankings.csv")))
pat=re.compile(r'extract|tea|juice|oil|powder|honey|jelly|ginseng|cranberry|apple|blueberry|grape|pomegranate|mushroom|herb|bark|root|leaf|seed|flower|fruit|berry|wine|propolis|chrysanthemum|astragalus|rhodiola|ashwagandha',re.I)
pure=[r for r in rankings if r["evidence_tier"]=="robust" and not pat.search(r["compound_name"])]
ext=[r for r in rankings if r["evidence_tier"]=="robust" and pat.search(r["compound_name"])]
out=ROOT/"outputs"/"canonical"/"compound_class_analysis.json"
out.write_text(json.dumps({"pure_robust":len(pure),"extract_robust":len(ext),"total_robust":len(pure)+len(ext)},indent=2)+"\n")
print(f"Wrote {out}")
