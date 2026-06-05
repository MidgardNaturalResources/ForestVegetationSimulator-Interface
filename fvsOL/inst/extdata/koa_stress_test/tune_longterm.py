"""tune_longterm.py -- calibrate the self-thinning ramp on 150-yr projections
and confirm long-term behavior for natural vs plantation across BYI."""
import os, numpy as np, pandas as pd
from koa_projector import project_cohort, SDI_MAX
from koa_survival_calibrated_py import make

OUT="results"; os.makedirs(OUT,exist_ok=True)
BYIS=[100,264,450]

def trajectory(byi, planted, ramp, max_age=150):
    fn = make(ramp=ramp)
    return project_cohort(byi, planted, "A", bounded=True, surv_fn=fn, max_age=max_age)

# --- sweep ramp coefficient: want SDI to peak near SDImax then track it, not blow past or collapse
print("="*78); print("RAMP SWEEP (natural BYI=264): peak SDI / age-150 SDI / age-150 QMD / age-150 TPH")
for ramp in (0.25,0.40,0.55,0.75,1.0):
    df=trajectory(264,0,ramp)
    print(f"  ramp={ramp:.2f}: peakSDI {df.SDI.max():4.0f} ({df.SDI.max()/SDI_MAX*100:3.0f}%) | "
          f"SDI150 {df.SDI.iloc[-1]:4.0f} | QMD150 {df.QMD.iloc[-1]:5.1f} | TPH150 {df.TPH.iloc[-1]:5.0f}")

RAMP=0.55   # selected below
print(f"\nSelected ramp={RAMP}")
print("="*78); print("LONG-TERM (150 yr) trajectories, ramp=%.2f" % RAMP)
rows=[]
for planted in (0,1):
    for byi in BYIS:
        df=trajectory(byi,planted,RAMP)
        for a in (20,40,60,100,150):
            r=df[df.age==a].iloc[0]
            rows.append(dict(origin="plt" if planted else "nat",BYI=byi,age=a,
                QMD=round(r.QMD,1),HT=round(r.HT,1),BAPH=round(r.BAPH,1),
                TPH=round(r.TPH,0),SDI=round(r.SDI,0),pctSDImax=round(r.SDI/SDI_MAX*100,0)))
lt=pd.DataFrame(rows); lt.to_csv(f"{OUT}/longterm_calibrated.csv",index=False)
for planted in (0,1):
    o="plt" if planted else "nat"
    print(f"\n {'PLANTATION' if planted else 'NATURAL'}:")
    sub=lt[lt.origin==o]
    for byi in BYIS:
        s=sub[sub.BYI==byi]
        traj=" ".join(f"a{r.age}:QMD{r.QMD:.0f}/TPH{r.TPH:.0f}/SDI{r.SDI:.0f}" for _,r in s.iterrows())
        print(f"  BYI{byi}: {traj}")

# --- check: does BYI act through growth (higher BYI reaches self-thinning sooner)? ---
print("\n"+"="*78); print("BYI-through-growth check: age stand first exceeds 55% SDImax (self-thin onset)")
for planted in (0,1):
    for byi in BYIS:
        df=trajectory(byi,planted,RAMP)
        on=df[df.SDI>=0.55*SDI_MAX]
        age_on=int(on.age.iloc[0]) if len(on) else None
        print(f"  {'plt' if planted else 'nat'} BYI{byi}: onset age {age_on}, peak QMD {df.QMD.max():.0f} cm")

# --- annual mortality realism at low density (background) ---
from koa_survival_calibrated_py import surv_calibrated
print("\n"+"="*78); print("Background annual mortality (RD=0.2): nat",
      round(1-surv_calibrated(15,10,0.5,0.5,264,sdi=100,planted=0),4),
      " plt", round(1-surv_calibrated(15,10,0.5,0.5,264,sdi=100,planted=1),4),
      "| near SDImax (RD=1.0): ", round(1-surv_calibrated(30,20,0.4,0.5,264,sdi=500,planted=0),4))
print("\nWrote longterm_calibrated.csv")
