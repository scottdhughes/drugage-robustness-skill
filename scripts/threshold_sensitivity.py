#!/usr/bin/env python3
"""Threshold sensitivity sweep."""
import csv,json,numpy as np
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
profiles=list(csv.DictReader(open(ROOT/"outputs"/"canonical"/"compound_evidence_profiles.csv")))
data=[{"name":p["compound_name"].strip(),"n_exp":int(p["num_experiments"]),"n_sp":int(p["num_species"]),
    "n_pmids":int(p["num_pmids"]),"sc":float(p["sign_consistency"]),"loso":float(p["leave_one_species_out_stability"]),
    "ag":float(p["aggregation_stability"]),"med":float(p["median_effect"]),"tm":float(p["trimmed_mean_effect"])} for p in profiles]
canon=set(d["name"] for d in data if d["n_exp"]>=3 and d["n_sp"]>=2 and d["n_pmids"]>=2 and d["sc"]>=0.80 and d["loso"]==1.0 and d["ag"]==1.0 and d["med"]>0 and d["tm"]>0)
results=[]
for me in [2,3,4,5]:
 for ms in [1,2,3]:
  for mp in [1,2,3]:
   for msc in [0.70,0.75,0.80,0.85,0.90]:
    for lr in [0.8,0.9,1.0]:
     rob=set(d["name"] for d in data if d["n_exp"]>=me and d["n_sp"]>=ms and d["n_pmids"]>=mp and d["sc"]>=msc and d["loso"]>=lr and d["ag"]==1.0 and d["med"]>0 and d["tm"]>0)
     ov=len(rob&canon); j=ov/max(len(rob|canon),1)
     results.append({"j":round(j,4),"n":len(rob)})
js=[r["j"] for r in results]
out=ROOT/"outputs"/"canonical"/"threshold_sensitivity.json"
out.write_text(json.dumps({"n_combinations":len(results),"tau_median":round(float(np.median(js)),4),
    "tau_mean":round(float(np.mean(js)),4),"fraction_above_050":round(sum(1 for j in js if j>0.50)/len(js),4),
    "fraction_above_070":round(sum(1 for j in js if j>0.70)/len(js),4),
    "robust_count_range":[min(r["n"] for r in results),max(r["n"] for r in results)]},indent=2)+"\n")
print(f"Wrote {out}")
