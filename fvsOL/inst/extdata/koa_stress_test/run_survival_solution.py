"""run_survival_solution.py
Select a survival solution by PROJECTION behavior (not just AUC), with focus on
plantations. Compares: S1 raw (current HiGy.R), S1 stabilized (clamp+floor),
and refit candidates, on (a) CFI first->last validation and (b) synthetic
plantation/natural cohorts.
"""
import os, json, numpy as np, pandas as pd
from koa_equations import predict_HT, LineageA as A
from koa_projector import project_psp, sdi_of
from survival_candidates import make_surv_fn, predict as scpred, VARIANTS

DATA="/sessions/elegant-lucid-cannon/mnt/Koa"; OUT="results"; os.makedirs(OUT,exist_ok=True)
BYI=264

# candidate survival functions (name -> surv_fn(dbh,ht,cr,rht,byi,bal,baph,planted))
def s1_raw(dbh,ht,cr,rht,byi,bal=0,baph=0,planted=0,**k):
    return A.surv_annual(dbh,ht,cr,rht,byi)
def s1_stable(dbh,ht,cr,rht,byi,bal=0,baph=0,planted=0,**k):
    return A.surv_annual_stable(dbh,ht,cr,rht,byi)
def s_calibrated(dbh,ht,cr,rht,byi,bal=0,baph=0,planted=0,**k):
    # Calibrated mortality: anchor to OBSERVED annual rates by origin, modulate by
    # relative density (self-thinning onset ~55% of BA envelope) and small size.
    dbh=np.asarray(dbh,float)
    base = 0.005 if planted else 0.015          # observed annual mortality
    RD = float(baph)/60.0                        # BA envelope ~60 m2/ha (FIA Q90+)
    comp = 1.0 + 2.5*max(0.0, RD-0.55)/0.45      # competition lift past self-thin
    size = 1.0 + 0.8*np.exp(-dbh/8.0)            # small-tree lift
    mort = np.clip(base*comp*size, 0.002, 0.12)
    return 1.0 - mort
CANDS={"S1_raw(current)":s1_raw, "S1_stabilized":s1_stable,
       "Snew_logit(refit)":make_surv_fn("Snew_logit"),
       "S_calibrated":s_calibrated}

# --- headline: annual mortality for a typical small plantation vs natural tree
print("="*86)
print("Annual mortality for a typical small tree (DBH8 HT6 CR0.6 rHT0.5 BAL5 BAPH15 BYI264)")
print("  observed annual mortality ~ nat 1.6%, plt 0.5%")
for nm,fn in CANDS.items():
    mp=1-float(np.atleast_1d(fn(8.,6.,0.6,0.5,BYI,bal=5,baph=15,planted=1))[0])
    mn=1-float(np.atleast_1d(fn(8.,6.,0.6,0.5,BYI,bal=5,baph=15,planted=0))[0])
    print(f"  {nm:18s} natural {mn*100:5.1f}%/yr   plantation {mp*100:5.1f}%/yr")

# --- CFI first->last validation ---
psp=pd.read_csv(f"{DATA}/AK_PSP_data_compiled.csv"); psp=psp[(psp.DBH_cm>0)&psp.TPA.notna()]
plots=[]
for pid,gp in psp.groupby("PSP"):
    seqs=sorted(gp.MeasSeq.unique())
    if len(seqs)<2: continue
    g0=gp[gp.MeasSeq==seqs[0]]; g1=gp[gp.MeasSeq==seqs[-1]]
    yrs=g1.Yfloor.max()-g0.Yfloor.max()
    if pd.isna(yrs) or yrs<3: continue
    yrs=int(round(yrs))
    b0=(g0.DBH_cm**2*0.00007854*g0.TPA).sum(); t0=g0.TPA.sum(); q0=np.sqrt(b0/(0.00007854*t0)) if t0>0 else 1
    tl=pd.DataFrame(dict(dbh=g0.DBH_cm.values,ht=np.where(g0.HT_m.values>1.37,g0.HT_m.values,np.nan),expf=g0.TPA.values))
    tl["ht"]=tl.ht.fillna(pd.Series(predict_HT(tl.dbh.values,max(b0,0.1),max(q0,1),BYI))); tl["cr"]=0.5
    bo=(g1.DBH_cm**2*0.00007854*g1.TPA).sum(); to=g1.TPA.sum(); qo=np.sqrt(bo/(0.00007854*to)) if to>0 else np.nan
    plots.append((pid,yrs,tl,qo,bo,to))
print("\n"+"="*86); print(f"CFI first->last validation (n={len(plots)} plots): bias vs observed")
cfi_rows=[]
for nm,fn in CANDS.items():
    QP,BP,TP,QO,BO,TO=[],[],[],[],[],[]
    for pid,yrs,tl,qo,bo,to in plots:
        p=project_psp(tl,byi=BYI,planted=0,lineage="A",n_years=yrs,surv_fn=fn)
        QP.append(p["QMD"]);BP.append(p["BAPH"]);TP.append(p["TPH"]);QO.append(qo);BO.append(bo);TO.append(to)
    QP,BP,TP,QO,BO,TO=map(np.array,(QP,BP,TP,QO,BO,TO))
    row=dict(candidate=nm,
        QMD_bias_pct=round(100*(QP-QO).mean()/QO.mean(),1),
        BAPH_bias_pct=round(100*(BP-BO).mean()/BO.mean(),1),
        TPH_bias_pct=round(100*(TP-TO).mean()/TO.mean(),1))
    cfi_rows.append(row)
    print(f"  {nm:18s} QMD {row['QMD_bias_pct']:6.1f}%  BAPH {row['BAPH_bias_pct']:7.1f}%  TPH {row['TPH_bias_pct']:7.1f}%")
pd.DataFrame(cfi_rows).to_csv(f"{OUT}/survival_solution_cfi.csv",index=False)

# --- synthetic cohorts: dense plantation & natural, 40 yr survival ---
def synth(planted):
    rng=np.random.default_rng(7)
    n=30
    dbh=np.clip(rng.normal(8 if planted else 12, 2.5, n),2,30)
    tph_total=1200 if planted else 500
    tl=pd.DataFrame(dict(dbh=dbh,
        ht=predict_HT(dbh,15.0,np.sqrt((dbh**2).mean()),BYI),
        expf=np.full(n,tph_total/n), cr=0.6))
    return tl
print("\n"+"="*86); print("Cohort 40-yr survival (% of initial TPH remaining); observed implies high retention")
coh_rows=[]
for origin,planted in (("natural",0),("plantation",1)):
    tl=synth(planted); tph0=tl.expf.sum()
    line=[f"  {origin:11s}"]
    rec=dict(origin=origin)
    for nm,fn in CANDS.items():
        p=project_psp(tl.copy(),byi=BYI,planted=planted,lineage="A",n_years=40,surv_fn=fn)
        pct=100*p["TPH"]/tph0; rec[nm]=round(pct,1); line.append(f"{nm.split('(')[0]}:{pct:5.1f}%")
    coh_rows.append(rec); print("  ".join(line))
pd.DataFrame(coh_rows).to_csv(f"{OUT}/survival_solution_cohort.csv",index=False)

json.dump(dict(cfi=cfi_rows,cohort=coh_rows),open(f"{OUT}/survival_solution_summary.json","w"),indent=2,default=str)
print("\nWrote survival_solution_cfi.csv, survival_solution_cohort.csv, summary json")
