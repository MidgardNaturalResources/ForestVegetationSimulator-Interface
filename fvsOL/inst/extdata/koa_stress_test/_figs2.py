import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import pandas as pd
from koa_equations import LineageA as A
fig,ax=plt.subplots(1,2,figsize=(11,4.2))
# (1) survival surface vs CR for rHT
crs=np.linspace(0.2,0.85,60)
for rht,c in ((0.3,"#2166ac"),(0.5,"#1a9850"),(0.7,"#d73027")):
    ax[0].plot(crs,[float(A.surv_annual(12,8,cr,rht,264)) for cr in crs],color=c,lw=2,label=f"rHT={rht}")
ax[0].axvspan(0.20,0.55,color="grey",alpha=0.15,label="stable clamp")
ax[0].set_xlabel("Crown ratio"); ax[0].set_ylabel("Annual survival prob")
ax[0].set_title("HiGy.R survival vs CR (DBH12,HT8,BYI264)"); ax[0].legend(fontsize=8)
# (2) CFI TPH bias raw vs stable
cfi=pd.read_csv("results/cfi_validation.csv")
xb=np.arange(3); 
raw=[ (cfi[f'{v}_raw']-cfi[f'{v}o']).mean()/cfi[f'{v}o'].mean()*100 for v in ("QMD","BAPH","TPH")]
stb=[ (cfi[f'{v}_stable']-cfi[f'{v}o']).mean()/cfi[f'{v}o'].mean()*100 for v in ("QMD","BAPH","TPH")]
ax[1].bar(xb-0.18,raw,0.36,label="raw HiGy.R",color="#d73027")
ax[1].bar(xb+0.18,stb,0.36,label="stabilized",color="#1a9850")
ax[1].set_xticks(xb); ax[1].set_xticklabels(["QMD","BAPH","TPH"]); ax[1].axhline(0,color="k",lw=0.7)
ax[1].set_ylabel("Mean bias (%)"); ax[1].set_title(f"CFI first->last bias (n={len(cfi)})"); ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig("results/fig_survival_dialin.png",dpi=130); print("fig ok")
