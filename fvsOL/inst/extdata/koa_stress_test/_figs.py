import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from koa_projector import project_cohort
R="results"
fig, ax = plt.subplots(1,2, figsize=(11,4.2))
for lin,col in (("A","#d73027"),("B","#2166ac")):
    for bounded,ls in ((True,"-"),(False,"--")):
        df=project_cohort(264,0,lin,bounded=bounded,max_age=100)
        ax[0].plot(df.age, df.QMD, ls, color=col, lw=1.8,
                   label=f"Lin {lin} {'bounded' if bounded else 'raw/FVS-like'}")
ax[0].axhline(42, color="k", ls=":", lw=1); ax[0].text(42,43,"Table 8 age-40 ~42cm",fontsize=8)
ax[0].set_xlabel("Stand age (yr)"); ax[0].set_ylabel("QMD (cm)")
ax[0].set_title("Natural cohort, BYI=264: QMD trajectory"); ax[0].legend(fontsize=8)
# SDI envelope
for lin,col in (("A","#d73027"),("B","#2166ac")):
    df=project_cohort(264,0,lin,bounded=True,max_age=100)
    ax[1].plot(df.age, df.SDI, color=col, lw=1.8, label=f"Lin {lin} (bounded)")
ax[1].axhline(500,color="k",ls=":"); ax[1].text(5,510,"SDImax=500",fontsize=8)
ax[1].set_xlabel("Stand age (yr)"); ax[1].set_ylabel("SDI"); ax[1].set_title("SDI vs envelope")
ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig(f"{R}/fig_cohort_trajectories.png", dpi=130)
print("figure written")
