import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
m=pd.read_csv("results/survival_instability_map.csv")
piv=m.pivot(index="rHT",columns="CR",values="surv")
fig,ax=plt.subplots(figsize=(6.2,3.6))
im=ax.imshow(piv.values,origin="lower",aspect="auto",cmap="RdYlGn",vmin=0,vmax=1,
   extent=[piv.columns.min(),piv.columns.max(),piv.index.min(),piv.index.max()])
ax.set_xlabel("Crown ratio"); ax.set_ylabel("Relative height (ht/htmax)")
ax.set_title("HiGy.R annual survival (DBH15, HT9, BYI264)")
cb=fig.colorbar(im,ax=ax); cb.set_label("P(survive)/yr")
ax.axvline(0.55,color="k",ls="--",lw=1); ax.text(0.56,0.32,"clamp 0.55",fontsize=8)
plt.tight_layout(); plt.savefig("results/fig_survival_heatmap.png",dpi=130); print("heatmap ok")
