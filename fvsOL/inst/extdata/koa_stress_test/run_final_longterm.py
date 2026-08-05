"""run_final_longterm.py -- FINAL long-term stress test of the calibrated survival
(linear self-thinning ramp 0.65->0.85 RD). Calibrate ramp height, then a full
nat/plt x BYI grid to 200 yr with robustness checks."""
import os, numpy as np, pandas as pd
from koa_projector import project_cohort, SDI_MAX
from koa_survival_calibrated_py import make, surv_calibrated
OUT="results"; os.makedirs(OUT,exist_ok=True); BYIS=[100,264,450]

def traj(byi,planted,maxlift,max_age=200):
    return project_cohort(byi,planted,"A",bounded=True,surv_fn=make(maxlift=maxlift),max_age=max_age)

print("="*80); print("RAMP-HEIGHT CALIBRATION (onset .65 -> full .85 RD), natural BYI=264")
print(" maxlift  peakRD  age@peak  RD200  QMD200  TPH200")
for ml in (0.04,0.06,0.08,0.10,0.12):
    df=traj(264,0,ml); pk=df.SDI.max()/SDI_MAX
    print(f"  {ml:.2f}    {pk:5.2f}    {int(df.loc[df.SDI.idxmax(),'age']):4d}     "
          f"{df.SDI.iloc[-1]/SDI_MAX:5.2f}  {df.QMD.iloc[-1]:5.1f}  {df.TPH.iloc[-1]:5.0f}")

ML=0.15
print(f"\nSelected maxlift={ML} (annual self-thinning mortality at >=85% RD)")
print("="*80); print("FINAL LONG-TERM GRID (200 yr): QMD / TPH / %SDImax")
rows=[]; checks=[]
for planted in (0,1):
    for byi in BYIS:
        df=traj(byi,planted,ML)
        peak=df.SDI.max()/SDI_MAX
        runaway = bool((df.SDI>1.05*SDI_MAX).any())
        collapse = bool((df.TPH<2).any())
        mono_qmd = bool((df.QMD.diff().dropna()>=-1e-6).all())  # QMD non-decreasing
        checks.append(dict(origin="plt" if planted else "nat",BYI=byi,
            peak_pctSDImax=round(peak*100,0), runaway=runaway, collapse=collapse,
            qmd_monotonic=mono_qmd, QMD200=round(df.QMD.iloc[-1],1), TPH200=round(df.TPH.iloc[-1],0)))
        for a in (25,50,100,150,200):
            r=df[df.age==a].iloc[0]
            rows.append(dict(origin="plt" if planted else "nat",BYI=byi,age=a,
                QMD=round(r.QMD,1),HT=round(r.HT,1),BAPH=round(r.BAPH,1),
                TPH=round(r.TPH,0),SDI=round(r.SDI,0),pctSDImax=round(r.SDI/SDI_MAX*100,0)))
lt=pd.DataFrame(rows); lt.to_csv(f"{OUT}/final_longterm.csv",index=False)
chk=pd.DataFrame(checks); chk.to_csv(f"{OUT}/final_longterm_checks.csv",index=False)
for planted in (0,1):
    o="plt" if planted else "nat"; print(f"\n {'PLANTATION' if planted else 'NATURAL'}:")
    for byi in BYIS:
        s=lt[(lt.origin==o)&(lt.BYI==byi)]
        print(f"  BYI{byi}: "+"  ".join(f"a{r.age}:{r.QMD:.0f}cm/{r.TPH:.0f}tph/{r.pctSDImax:.0f}%" for _,r in s.iterrows()))
print("\n"+"="*80); print("ROBUSTNESS CHECKS (all should be: no runaway, no collapse, QMD monotonic)")
print(chk.to_string(index=False))
print(f"\n  ANY runaway: {chk.runaway.any()} | ANY collapse: {chk.collapse.any()} | "
      f"ALL QMD monotonic: {chk.qmd_monotonic.all()} | peak %SDImax range "
      f"{chk.peak_pctSDImax.min():.0f}-{chk.peak_pctSDImax.max():.0f}")

# background + plateau annual mortality sanity
print("\nAnnual mortality: background(RD .2) nat",
      round(1-surv_calibrated(15,10,.5,.5,264,sdi=0.2*SDI_MAX,planted=0,maxlift=ML),4),
      "plt",round(1-surv_calibrated(15,10,.5,.5,264,sdi=0.2*SDI_MAX,planted=1,maxlift=ML),4),
      "| at full self-thin (RD .9) nat",
      round(1-surv_calibrated(30,20,.4,.5,264,sdi=0.9*SDI_MAX,planted=0,maxlift=ML),4))
print("Wrote final_longterm.csv, final_longterm_checks.csv")
