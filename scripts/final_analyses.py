#!/usr/bin/env python3
"""Unseen-species, pure-compound, and simple baseline analyses."""
import csv,json,re,numpy as np
from collections import defaultdict
from pathlib import Path
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score
ROOT=Path(__file__).resolve().parents[1]
rankings=list(csv.DictReader(open(ROOT/"outputs"/"canonical"/"robustness_rankings.csv")))
cr={r["compound_name"].strip():int(r["rank"]) for r in rankings}; tqc=len(rankings)//4; tq=set(c for c,r in cr.items() if r<=tqc)
itp=json.loads((ROOT/"outputs"/"canonical"/"itp_overlap_analysis.json").read_text())
ip=set(d["compound"] for d in itp["compound_details"] if d["itp_direction"]=="positive")
ine=set(d["compound"] for d in itp["compound_details"] if d["itp_direction"]=="negative")
hr=list(csv.DictReader(open(ROOT/"outputs"/"itp_holdout"/"run"/"robustness_rankings.csv")))
hm={r["compound_name"].strip():float(r["robustness_score"]) for r in hr}
hp=list(csv.DictReader(open(ROOT/"outputs"/"itp_holdout"/"run"/"compound_evidence_profiles.csv")))
hs={}
for p in hp:
    n=p["compound_name"].strip();ns=int(p["num_species"]);sc=float(p["sign_consistency"]);md=float(p["median_effect"])
    hs[n]=ns*sc if md>0 else 0
ae=sorted((ip|ine)&set(hm)&set(hs)); lb=[1 if c in ip else 0 for c in ae]
ar=roc_auc_score(lb,[hm[c] for c in ae]); asb=roc_auc_score(lb,[hs[c] for c in ae])
# Unseen species
b4s=defaultdict(set)
with open(ROOT/"data"/"drugage_build4_2021-11-20.csv",encoding="utf-8",errors="replace") as f:
    for r in csv.DictReader(f): b4s[r.get("compound_name","").strip()].add(r.get("species","").strip())
ut,ub=[],[]
b5c=defaultdict(list)
with open(ROOT/"data"/"drugage_build5_2024-11-29.csv",encoding="utf-8") as f:
    for r in csv.DictReader(f):
        c=r.get("compound_name","").strip();sp=r.get("species","").strip()
        try:
            e=float(r["avg_lifespan_change_percent"])
            if c in b4s and sp not in b4s[c] and c in cr: b5c[c].append(1 if e>0 else 0)
        except: pass
for c,o in b5c.items():
    r=np.mean(o)
    if c in tq: ut.append(r)
    else: ub.append(r)
up=mannwhitneyu(ut,ub,alternative="greater")[1] if ut and ub else 1
result={"simple_baseline_auroc":round(float(asb),4),"robustness_auroc":round(float(ar),4),
    "unseen_top_q_rate":round(float(np.mean(ut)),4) if ut else None,
    "unseen_bottom_rate":round(float(np.mean(ub)),4) if ub else None,"unseen_p":round(float(up),6)}
out=ROOT/"outputs"/"canonical"/"final_analyses.json"
out.write_text(json.dumps(result,indent=2)+"\n"); print(f"Wrote {out}")
