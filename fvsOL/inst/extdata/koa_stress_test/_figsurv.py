import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import pandas as pd
cfi=pd.read_csv("results/survival_solution_cfi.csv")
coh=pd.read_csv("results/survival_solution_cohort.csv")
cands=["S1_raw(current)","S1_stabilized","Snew_logit(refit)","S_calibrated"]
short=["S1 raw\n(current)","S1\nstabilized","refit\n(overfit)","calibrated\n(recommended)"]
fig,ax=plt.subplots(1,2,figsize=(11,4.3))
# panel1: 40-yr cohort retention
nat=[float(coh[coh.origin=="natural"][c]) for c in cands]
plt_=[float(coh[coh.origin=="plantation"][c]) for c in cands]
x=np.arange(len(cands))
ax[0].bar(x-0.18,nat,0.36,label="natural",color="#1a9850")
ax[0].bar(x+0.18,plt_,0.36,label="plantation",color="#d73027")
ax[0].axhspan(50,85,color="grey",alpha=0.15,label="realistic band")
ax[0].set_xticks(x); ax[0].set_xticklabels(short,fontsize=8); ax[0].set_ylabel("40-yr TPH retained (%)")
ax[0].set_title("Cohort survival behavior"); ax[0].legend(fontsize=8)
# panel2: typical-tree annual mortality
mort={"S1_raw(current)":(100,100),"S1_stabilized":(10,10),"Snew_logit(refit)":(0,99.2),"S_calibrated":(1.9,0.6)}
mn=[mort[c][0] for c in cands]; mp=[mort[c][1] for c in cands]
ax[1].bar(x-0.18,mn,0.36,label="natural",color="#1a9850")
ax[1].bar(x+0.18,mp,0.36,label="plantation",color="#d73027")
ax[1].axhline(1.6,color="#1a9850",ls=":",lw=1.2); ax[1].axhline(0.5,color="#d73027",ls=":",lw=1.2)
ax[1].set_xticks(x); ax[1].set_xticklabels(short,fontsize=8); ax[1].set_ylabel("Annual mortality, typical small tree (%)")
ax[1].set_title("Annual mortality vs observed (dotted)"); ax[1].set_ylim(0,105); ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig("results/fig_survival_solution.png",dpi=130); print("fig ok")
