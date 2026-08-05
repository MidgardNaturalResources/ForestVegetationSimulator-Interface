"""run_origin_byi_final.py -- (1) reconfirm long-term robustness with corrected
origin base rates; (2) show that 'survival is lower at higher BYI' EMERGES from
the growth -> self-thinning pathway (no direct BYI mortality term)."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from koa_projector import project_cohort, SDI_MAX
from koa_survival_calibrated_py import make
OUT="results"; fn=make(); BYIS=[100,150,200,264,350,450,550]

# robustness with corrected bases
print("Robustness (corrected bases nat .003/plt .006):")
for planted in (0,1):
    for byi in (100,264,450):
        df=project_cohort(byi,planted,"A",bounded=True,surv_fn=fn,max_age=200)
        print(f"  {'plt' if planted else 'nat'} BYI{byi}: peak {df.SDI.max()/SDI_MAX*100:.0f}%SDImax "
              f"QMD200 {df.QMD.iloc[-1]:.0f} TPH200 {df.TPH.iloc[-1]:.0f} "
              f"runaway={(df.SDI>1.05*SDI_MAX).any()} collapse={(df.TPH<2).any()}")

# emergent BYI -> mortality (natural): mean annual mortality & cumulative survival to 100
print("\nEmergent realized mortality by BYI (natural), age 1-100:")
rows=[]
for planted in (0,1):
    for byi in BYIS:
        df=project_cohort(byi,planted,"A",bounded=True,surv_fn=fn,init_dbh=2.0,
                          init_tph=(1500 if planted else 800),init_age=1,max_age=100)
        tph=df.TPH.values
        cum_surv=tph[-1]/tph[0]
        ann=df.SDI.values  # placeholder
        # realized mean annual mortality = 1 - (cum_surv)^(1/years)
        yrs=df.age.iloc[-1]-df.age.iloc[0]
        mean_ann=1-cum_surv**(1/max(yrs,1))
        onset=df[df.SDI>=0.65*SDI_MAX]
        rows.append(dict(origin="plt" if planted else "nat",BYI=byi,
            mean_ann_mort=round(mean_ann,4), cum_surv_100=round(cum_surv,3),
            selfthin_onset_age=int(onset.age.iloc[0]) if len(onset) else None,
            QMD100=round(df.QMD.iloc[-1],0)))
em=pd.DataFrame(rows); em.to_csv(f"{OUT}/emergent_byi_mortality.csv",index=False)
for o in ("nat","plt"):
    s=em[em.origin==o]
    print(f" {o}: "+"  ".join(f"BYI{r.BYI}:{r.mean_ann_mort*100:.1f}%/yr(onset~a{r.selfthin_onset_age})" for _,r in s.iterrows()))

# figure: realized mean annual mortality vs BYI (emergent), both origins
fig,ax=plt.subplots(1,2,figsize=(11,4.2))
for o,c in (("nat","#1a9850"),("plt","#d73027")):
    s=em[em.origin==o]
    ax[0].plot(s.BYI,s.mean_ann_mort*100,"o-",color=c,label=o)
    ax[1].plot(s.BYI,s.selfthin_onset_age,"o-",color=c,label=o)
ax[0].set_xlabel("BYI (Mg/ha)"); ax[0].set_ylabel("Realized mean annual mortality, age 1-100 (%)")
ax[0].set_title("Survival lower at higher BYI (emergent via self-thinning)"); ax[0].legend(fontsize=8)
ax[1].set_xlabel("BYI (Mg/ha)"); ax[1].set_ylabel("Self-thinning onset age (yr)")
ax[1].set_title("Higher BYI reaches self-thinning sooner"); ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig(f"{OUT}/fig_emergent_byi.png",dpi=130)
print("\nWrote emergent_byi_mortality.csv, fig_emergent_byi.png")
