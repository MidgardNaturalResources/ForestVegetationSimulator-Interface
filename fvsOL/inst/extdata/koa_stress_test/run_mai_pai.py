"""run_mai_pai.py -- PAI/MAI volume trends and MAI culmination (biological
rotation age) from the final calibrated long-term projections."""
import os, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from koa_projector import project_cohort
from koa_survival_calibrated_py import make
OUT="results"; BYIS=[100,264,450]; fn=make()

def series(byi,planted,max_age=200):
    # Start from establishment (small DBH, age 1) so MAI/PAI and culmination are
    # measured over true stand age from regeneration/planting.
    df=project_cohort(byi,planted,"A",bounded=True,surv_fn=fn,
                      init_dbh=2.0, init_tph=(1500 if planted else 800),
                      init_age=1, max_age=max_age)
    age=df.age.values.astype(float); vol=df.VOL.values
    mai=vol/age
    pai=np.gradient(vol, age)            # annual volume increment
    return age,vol,mai,pai

rows=[]; fig,ax=plt.subplots(1,2,figsize=(11,4.3))
col={100:"#2166ac",264:"#1a9850",450:"#d73027"}
for planted,ls,axi in ((0,"-",0),(1,"--",1)):
    for byi in BYIS:
        age,vol,mai,pai=series(byi,planted)
        cul=int(age[np.argmax(mai)]); maipk=float(np.max(mai))
        # MAI=PAI crossover (first age after peak where pai<=mai)
        cross=age[(age>cul)&(pai<=mai)]
        cross=int(cross[0]) if len(cross) else None
        rows.append(dict(origin="plt" if planted else "nat",BYI=byi,
            MAI_culm_age=cul, MAI_peak_m3ha_yr=round(maipk,2),
            PAI_peak_age=int(age[np.argmax(pai)]), vol200=round(vol[-1],0)))
        ax[axi].plot(age,mai,ls,color=col[byi],lw=1.8,label=f"MAI BYI{byi}")
        ax[axi].plot(age,pai,ls,color=col[byi],lw=1.0,alpha=0.6)
        ax[axi].axvline(cul,color=col[byi],ls=":",lw=0.8)
for a,t in zip(ax,("Natural","Plantation")):
    a.set_xlabel("Stand age (yr)"); a.set_ylabel("Volume increment (m3/ha/yr)")
    a.set_title(f"{t}: MAI (bold) & PAI (thin); dotted = MAI culmination"); a.legend(fontsize=7)
plt.tight_layout(); plt.savefig(f"{OUT}/fig_mai_pai.png",dpi=130)
mp=pd.DataFrame(rows); mp.to_csv(f"{OUT}/mai_pai.csv",index=False)
print("MAI/PAI summary (volume = BAPH*HT*0.40):"); print(mp.to_string(index=False))
print("\nWrote mai_pai.csv, fig_mai_pai.png")
