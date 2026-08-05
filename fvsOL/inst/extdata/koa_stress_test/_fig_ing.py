import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from koa_ingrowth import ingrowth_annual
lt=pd.read_csv("results/ingrowth_longterm.csv")
fig,ax=plt.subplots(1,2,figsize=(11,4.2))
rd=np.linspace(0,1,50)
ax[0].plot(rd,[ingrowth_annual(sdi=r*500,planted=0) for r in rd],color="#1a9850",lw=2,label="natural")
ax[0].plot(rd,[ingrowth_annual(sdi=r*500,planted=1) for r in rd],color="#d73027",lw=2,label="plantation")
ax[0].set_xlabel("Relative density (SDI/SDImax)"); ax[0].set_ylabel("Ingrowth (trees/ha/yr)")
ax[0].set_title("Koa ingrowth = exp(5.38 - 3.09*RD - 1.64*planted)"); ax[0].legend(fontsize=8)
# long-term BA with/without ingrowth (natural)
for ig,c,lab in ((False,"#888888","no ingrowth"),(True,"#1a9850","with ingrowth")):
    s=lt[(lt.origin=="nat")&(lt.ingrowth==ig)].sort_values("age")
    ax[1].plot(s.age,s.BAPH,"o-",color=c,lw=1.8,label=lab)
ax[1].set_xlabel("Stand age (yr)"); ax[1].set_ylabel("Basal area (m2/ha), natural")
ax[1].set_title("Ingrowth sustains long-term basal area"); ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig("results/fig_ingrowth.png",dpi=130); print("fig ok")
