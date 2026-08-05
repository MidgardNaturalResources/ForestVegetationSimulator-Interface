import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from koa_projector import project_cohort, SDI_MAX
from koa_survival_calibrated_py import make
fn=make()  # final defaults: onset .65, full .85, maxlift .15
fig,ax=plt.subplots(1,3,figsize=(13,4))
col={100:"#2166ac",264:"#1a9850",450:"#d73027"}
for planted,ls in ((0,"-"),(1,"--")):
    for byi in (100,264,450):
        df=project_cohort(byi,planted,"A",bounded=True,surv_fn=fn,max_age=200)
        lab=f"BYI{byi} {'plt' if planted else 'nat'}"
        ax[0].plot(df.age,df.SDI,ls,color=col[byi],lw=1.6,label=lab)
        ax[1].plot(df.age,df.QMD,ls,color=col[byi],lw=1.6)
        ax[2].plot(df.QMD,df.TPH,ls,color=col[byi],lw=1.6)
ax[0].axhline(SDI_MAX,color="k",ls=":"); ax[0].text(5,SDI_MAX+6,"SDImax",fontsize=8)
ax[0].axhline(0.65*SDI_MAX,color="grey",ls=":",lw=0.8); ax[0].text(5,0.65*SDI_MAX+6,"0.65 RD onset",fontsize=7)
ax[0].set_xlabel("Age (yr)"); ax[0].set_ylabel("SDI"); ax[0].set_title("SDI: self-thinning (ramp .65->.85 RD)")
ax[0].legend(fontsize=6,ncol=2)
ax[1].set_xlabel("Age (yr)"); ax[1].set_ylabel("QMD (cm)"); ax[1].set_title("QMD (solid=nat, dash=plt)")
ax[2].set_xlabel("QMD (cm)"); ax[2].set_ylabel("TPH"); ax[2].set_yscale("log"); ax[2].set_title("Stand trajectory (Bakuzis)")
plt.tight_layout(); plt.savefig("results/fig_longterm_calibrated.png",dpi=130); print("fig ok")
