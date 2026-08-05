"""run_ingrowth_behavior.py -- how ingrowth behaves over long-term (300 yr):
the annual recruitment rate trajectory and whether recruitment settles into
turnover balance with mortality at the steady state."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from koa_projector import project_psp, SDI_MAX
from koa_equations import predict_HT
np.random.seed(11)

def cohort(planted):
    dbh=np.clip(np.random.normal(8 if planted else 12,3,30),2,30); tph=1200 if planted else 600
    return pd.DataFrame(dict(dbh=dbh,ht=predict_HT(dbh,15.,np.sqrt((dbh**2).mean()),264),
                             expf=np.full(len(dbh),tph/len(dbh)),cr=0.6))

print("="*82); print("INGROWTH over 300 yr (trees/ha/yr): early -> decline -> steady turnover")
fig,ax=plt.subplots(1,2,figsize=(11,4.2)); col={100:"#2166ac",264:"#1a9850",450:"#d73027"}
rows=[]
for planted,ls in ((0,"-"),(1,"--")):
    for byi in (100,264,450):
        o=project_psp(cohort(planted),byi,planted,"A",300,surv_mode="calib_alloc",
                      ingrowth=True,return_traj=True); tr=o["traj"]; age=np.arange(1,len(tr)+1)
        # mortality rate (trees/ha/yr) = -d(TPH)/dt + ingrowth (since dTPH = ingrowth - mort)
        dT=np.diff(tr.TPH.values,prepend=tr.TPH.values[0])
        mort=tr.ingrowth.values - dT      # mort = ingrowth - change in TPH
        ax[0].plot(age,tr.ingrowth,ls,color=col[byi],lw=1.5,label=f"BYI{byi} {'plt' if planted else 'nat'}")
        early=tr.ingrowth.iloc[:10].mean(); late=tr.ingrowth.iloc[-50:].mean()
        late_mort=mort[-50:].mean()
        rows.append(dict(origin="plt" if planted else "nat",BYI=byi,
            ing_age1_10=round(early,1), ing_steady=round(late,1),
            mort_steady=round(late_mort,1), balance=round(late-late_mort,1),
            SDI300=round(tr.SDI.iloc[-1]/SDI_MAX*100)))
        if planted==0: ax[1].plot(age,tr.ingrowth,color=col[byi],lw=1.3,label=f"ingrowth BYI{byi}")
        if planted==0: ax[1].plot(age,mort,color=col[byi],lw=1.0,ls=":")
ax[0].set_xlabel("Age (yr)"); ax[0].set_ylabel("Ingrowth (trees/ha/yr)")
ax[0].set_title("Ingrowth trajectory (solid=nat, dash=plt)"); ax[0].legend(fontsize=6,ncol=2)
ax[1].set_xlabel("Age (yr)"); ax[1].set_ylabel("trees/ha/yr (natural)")
ax[1].set_title("Ingrowth (solid) vs mortality (dotted): turnover balance"); ax[1].legend(fontsize=7)
plt.tight_layout(); plt.savefig("results/fig_ingrowth_behavior.png",dpi=130)
b=pd.DataFrame(rows); b.to_csv("results/ingrowth_behavior.csv",index=False)
print(b.to_string(index=False))
print("\nInterpretation: ingrowth starts high in open stands, declines as RD rises,")
print("then settles to a steady rate that ~balances mortality (turnover) at equilibrium.")
print("Wrote ingrowth_behavior.csv, fig_ingrowth_behavior.png")
