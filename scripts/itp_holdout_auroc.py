#!/usr/bin/env python3
"""ITP holdout AUROC with bootstrap CIs."""
import csv,json,numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score
ROOT=Path(__file__).resolve().parents[1]
hr=list(csv.DictReader(open(ROOT/"outputs"/"itp_holdout"/"run"/"robustness_rankings.csv")))
hm={r["compound_name"].strip():float(r["robustness_score"]) for r in hr}
itp=json.loads((ROOT/"outputs"/"canonical"/"itp_overlap_analysis.json").read_text())
pos=set(d["compound"] for d in itp["compound_details"] if d["itp_direction"]=="positive")
neg=set(d["compound"] for d in itp["compound_details"] if d["itp_direction"]=="negative")
ae=sorted((pos|neg)&set(hm)); labels=np.array([1 if c in pos else 0 for c in ae])
scores=np.array([hm[c] for c in ae]); auc=roc_auc_score(labels,scores)
rng=np.random.default_rng(42); aucs=[]
for _ in range(2000):
    idx=rng.choice(len(labels),size=len(labels),replace=True)
    if len(set(labels[idx]))<2: continue
    aucs.append(roc_auc_score(labels[idx],scores[idx]))
aucs=np.array(aucs)
result={"n_positive":int(sum(labels)),"n_negative":int(len(labels)-sum(labels)),
    "auroc_by_metric":{"robustness_score":round(float(auc),4)},
    "auroc_bootstrap_ci_low":round(float(np.percentile(aucs,2.5)),4),
    "auroc_bootstrap_ci_high":round(float(np.percentile(aucs,97.5)),4),
    "auroc_bootstrap_n":len(aucs)}
out=ROOT/"outputs"/"canonical"/"itp_holdout_auroc.json"
out.write_text(json.dumps(result,indent=2)+"\n"); print(f"Wrote {out}")
