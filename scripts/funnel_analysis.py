#!/usr/bin/env python3
"""Positivity-skew audit."""
import csv,json,numpy as np
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
effects=[]
with open(ROOT/"data"/"drugage_build5_2024-11-29.csv",encoding="utf-8") as f:
    for r in csv.DictReader(f):
        try: effects.append(float(r["avg_lifespan_change_percent"]))
        except: pass
np_=sum(1 for e in effects if e>0); nn=len(effects)-np_
result={"total_experiments":len(effects),"positive_experiments":np_,"negative_experiments":nn,
    "positive_fraction":round(np_/len(effects),4),"median_effect":round(float(np.median(effects)),2)}
out=ROOT/"outputs"/"canonical"/"funnel_analysis.json"
out.write_text(json.dumps(result,indent=2)+"\n"); print(f"Wrote {out}")
