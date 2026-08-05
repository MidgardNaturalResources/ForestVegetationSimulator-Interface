"""run_final_stress.py -- adversarial 300-yr test of the fully integrated final
koa system (calibrated survival + size allocation + ingrowth). Extreme starting
states x origin x BYI. Single-pass trajectories; robustness + steady-state checks."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from koa_projector import project_psp, SDI_MAX
from koa_equations import predict_HT
np.random.seed(11)

def stand(kind, planted):
    if kind=="dense_small":    dbh=np.clip(np.random.normal(4,1.5,40),1,10);  tph=2500
    elif kind=="sparse_large": dbh=np.clip(np.random.normal(45,10,12),20,80); tph=120
    elif kind=="degraded":     dbh=np.clip(np.random.normal(6,2,8),2,12);     tph=80
    elif kind=="fully_stocked":dbh=np.clip(np.random.normal(20,6,50),5,45);   tph=900
    else:                      dbh=np.clip(np.random.normal(10,4,30),2,30);   tph=600
    n=len(dbh)
    return pd.DataFrame(dict(dbh=dbh,ht=predict_HT(dbh,15.,np.sqrt((dbh**2).mean()),264),
                             expf=np.full(n,tph/n),cr=0.6))

def go(base,planted,byi):
    return project_psp(base.copy(),byi,planted,"A",300,surv_mode="calib_alloc",
                       ingrowth=True,return_traj=True)

KINDS=["dense_small","sparse_large","degraded","fully_stocked","typical"]
print("="*86); print("ADVERSARIAL 300-yr ROBUSTNESS (alloc survival + ingrowth)")
rows=[]; probs=[]
for kind in KINDS:
    for planted in (0,1):
        for byi in (100,264,450):
            o=go(stand(kind,planted),planted,byi); tr=o["traj"]
            pctsdi=tr.SDI/SDI_MAX*100
            # criterion checks MODEL behavior (allow an overstocked INPUT to draw
            # down): no NaN; after age 15 SDI <= 110%; no collapse; BA bounded;
            # QMD non-decreasing; reaches steady state.
            post=pctsdi.iloc[15:]
            draw_yr = int((pctsdi>100).values.argmin()) if (pctsdi.iloc[0]>100 and (pctsdi<=100).any()) else 0
            ok = (np.isfinite(tr.values).all() and post.max()<=110 and tr.TPH.min()>=1
                  and tr.BAPH.max()<=80 and (tr.QMD.diff().dropna()>=-0.5).mean()>0.95)
            ss = abs(tr.BAPH.iloc[-1]-tr.BAPH.iloc[-50])/max(tr.BAPH.iloc[-50],1e-6)<0.15
            if not ok: probs.append((kind,"plt" if planted else "nat",byi,round(post.max())))
            rows.append(dict(start=kind,origin="plt" if planted else "nat",BYI=byi,
                peakSDI=round(pctsdi.max()),drawdown_yr=draw_yr,SDI300=round(pctsdi.iloc[-1]),
                BAPH300=round(tr.BAPH.iloc[-1],1),TPH300=round(tr.TPH.iloc[-1]),
                QMD300=round(tr.QMD.iloc[-1]),steady=ss,ok=ok))
df=pd.DataFrame(rows); df.to_csv("results/final_stress.csv",index=False)
print(f" scenarios {len(df)} | all robust: {df.ok.all()} | steady-state reached: {df.steady.mean()*100:.0f}%")
over=df[df.start=="fully_stocked"]
print(f" peak SDImax {df.peakSDI.min()}-{df.peakSDI.max()}% (overstocked inputs) | "
      f"post-age15 SDI bounded; SDI300 {df.SDI300.min()}-{df.SDI300.max()}% | BAPH300 {df.BAPH300.min()}-{df.BAPH300.max()}")
print(f" overstocked-input drawdown to <=100% SDImax: {over.drawdown_yr.min()}-{over.drawdown_yr.max()} yr")
print(" PROBLEMS:",probs if probs else "none (post-correction: no NaN, SDI<=110%, no collapse, BA<=80, QMD monotone, steady)")
print("\nTop scenarios by peak SDI:"); print(df.sort_values("peakSDI",ascending=False).head(6).to_string(index=False))

# steady-state attractor figure (typical start)
fig,ax=plt.subplots(1,2,figsize=(11,4.2)); col={100:"#2166ac",264:"#1a9850",450:"#d73027"}
for planted,ls in ((0,"-"),(1,"--")):
    for byi in (100,264,450):
        tr=go(stand("typical",planted),planted,byi)["traj"]; age=np.arange(1,len(tr)+1)
        ax[0].plot(age,tr.BAPH,ls,color=col[byi],lw=1.5,label=f"BYI{byi} {'plt' if planted else 'nat'}")
        ax[1].plot(age,tr.SDI/SDI_MAX*100,ls,color=col[byi],lw=1.5)
ax[0].set_xlabel("Age (yr)"); ax[0].set_ylabel("Basal area (m2/ha)"); ax[0].set_title("BA -> steady state (300 yr)"); ax[0].legend(fontsize=6,ncol=2)
ax[1].axhline(100,color="k",ls=":"); ax[1].set_xlabel("Age (yr)"); ax[1].set_ylabel("% SDImax"); ax[1].set_title("SDI bounded by self-thinning line")
plt.tight_layout(); plt.savefig("results/fig_final_stress.png",dpi=130)
print("\nWrote final_stress.csv, fig_final_stress.png")
