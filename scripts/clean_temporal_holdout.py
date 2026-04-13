#!/usr/bin/env python3
"""Clean Build-4-only temporal holdout. Zero leakage."""
import csv, json, numpy as np
from collections import defaultdict
from pathlib import Path
from scipy.stats import mannwhitneyu
import statsmodels.api as sm

ROOT = Path(__file__).resolve().parents[1]
WEIGHTS = np.array([0.35, 0.20, 0.15, 0.10, 0.15, 0.05])

def score(experiments):
    by_c = defaultdict(list)
    for e in experiments: by_c[e["compound"]].append(e)
    scores = {}
    for c, exps in by_c.items():
        n = len(exps); species = set(e["species"] for e in exps)
        n_pos = sum(e["positive"] for e in exps); effects = [e["effect"] for e in exps]
        sc = n_pos/n if n else 0; med = float(np.median(effects))
        tf = max(1,int(0.1*n)); se = sorted(effects)
        tm = float(np.mean(se[tf:n-tf])) if n>2 else float(np.mean(se))
        ag = 1.0 if med>0 and tm>0 else 0.0; mg = min(max(tm,0),50)/50
        if len(species)>=2:
            lp=0
            for sp in species:
                sub=[e["effect"] for e in exps if e["species"]!=sp]
                if sub and np.mean(sub)>0: lp+=1
            lo=lp/len(species)
        else: lo=0.0
        br=(min(n,6)/6+min(len(species),4)/4)/2
        scores[c]=float(np.array([br,sc,lo,0,ag,mg])@WEIGHTS)
    return scores

b4=[]
with open(ROOT/"data"/"drugage_build4_2021-11-20.csv",encoding="utf-8",errors="replace") as f:
    for r in csv.DictReader(f):
        try:
            e=float(r.get("avg_lifespan_change","0"))
            b4.append({"compound":r["compound_name"].strip(),"species":r["species"].strip(),"effect":e,"positive":1 if e>0 else 0})
        except: pass

b4s=score(b4); b4r=sorted(b4s,key=lambda c:-b4s[c]); tqc=len(b4r)//4; tq=set(b4r[:tqc])

b4c=defaultdict(list)
for e in b4: b4c[e["compound"]].append(e)
b5c=defaultdict(list)
with open(ROOT/"data"/"drugage_build5_2024-11-29.csv",encoding="utf-8") as f:
    for r in csv.DictReader(f):
        try:
            e=float(r["avg_lifespan_change_percent"])
            b5c[r["compound_name"].strip()].append({"effect":e,"positive":1 if e>0 else 0})
        except: pass

dd=[]
for c in b5c:
    b4n=len(b4c.get(c,[])); b5n=len(b5c[c])
    if b5n>b4n and c in b4s:
        de=b5c[c][b4n:]
        dd.append({"compound":c,"in_top_q":c in tq,"b4_n":b4n,"delta_n":b5n-b4n,"delta_pos_rate":np.mean([e["positive"] for e in de])})

tr=[d["delta_pos_rate"] for d in dd if d["in_top_q"]]
br2=[d["delta_pos_rate"] for d in dd if not d["in_top_q"]]
u,p=mannwhitneyu(tr,br2,alternative="greater")
Y=np.array([d["delta_pos_rate"] for d in dd])
X=sm.add_constant(np.column_stack([[d["in_top_q"] for d in dd],[np.log1p(d["b4_n"]) for d in dd]]))
m=sm.OLS(Y,X).fit()

# Propensity matching
top_d=[d for d in dd if d["in_top_q"]]; bot_d=[d for d in dd if not d["in_top_q"]]
mt,mb=[],[];used=set()
for t in sorted(top_d,key=lambda x:x["b4_n"]):
    best,bd=None,float("inf")
    for i,b in enumerate(bot_d):
        if i in used: continue
        diff=abs(t["b4_n"]-b["b4_n"])
        if diff<bd: bd=diff;best=i
    if best is not None and bd<=3: mt.append(t);mb.append(bot_d[best]);used.add(best)
mtr=[d["delta_pos_rate"] for d in mt]; mbr=[d["delta_pos_rate"] for d in mb]
u2,p2=mannwhitneyu(mtr,mbr,alternative="greater") if mt else (0,1)

result={"design":"Build-4-only scoring","b4_compounds":len(b4s),"b4_top_quartile":tqc,
    "compounds_with_delta":len(dd),"top_q_count":len(tr),"bottom_count":len(br2),
    "top_q_delta_pos_rate":round(float(np.mean(tr)),4),"bottom_delta_pos_rate":round(float(np.mean(br2)),4),
    "mann_whitney_p":round(float(p),6),"ols_top_q_coef":round(float(m.params[1]),4),
    "ols_top_q_p":round(float(m.pvalues[1]),6),"matched_pairs":len(mt),
    "matched_top_rate":round(float(np.mean(mtr)),4) if mt else None,
    "matched_bot_rate":round(float(np.mean(mbr)),4) if mb else None,
    "matched_p":round(float(p2),6) if mt else None}
out=ROOT/"outputs"/"canonical"/"clean_temporal_holdout.json"
out.write_text(json.dumps(result,indent=2)+"\n"); print(f"Wrote {out}")
