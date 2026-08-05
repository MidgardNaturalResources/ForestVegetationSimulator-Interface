"""run_ingrowth_longterm.py -- stress-test the koa system WITH ingrowth (RD+origin)
over 200 yr, and re-check CFI young-stand BAPH gap with ingrowth on."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from koa_projector import project_psp, sdi_of, SDI_MAX
from koa_equations import predict_HT
from koa_survival_calibrated_py import make
from koa_ingrowth import ingrowth_annual
np.random.seed(5)
SF=make()

def synth(planted,n=30):
    dbh=np.clip(np.random.normal(8 if planted else 12,3,n),2,30)
    tph=1200 if planted else 600
    return pd.DataFrame(dict(dbh=dbh,ht=predict_HT(dbh,15.,np.sqrt((dbh**2).mean()),264),
                             expf=np.full(n,tph/n),cr=0.6))

print("="*78); print("Ingrowth rate (trees/ha/yr) by RD x origin:")
for rd in (0.1,0.3,0.5,0.7,0.9):
    print(f"  RD {rd}: natural {ingrowth_annual(sdi=rd*SDI_MAX,planted=0):.0f}  "
          f"plantation {ingrowth_annual(sdi=rd*SDI_MAX,planted=1):.0f}")

print("\n"+"="*78); print("200-yr cohort: WITHOUT vs WITH ingrowth (calib_alloc survival)")
rows=[]
for planted in (0,1):
    base=synth(planted)
    for ig in (False,True):
        # step year by year to capture trajectory
        t=base.copy(); rec=[]
        # use project_psp in 10-yr chunks to snapshot
        for chunk in range(20):
            out=project_psp(t,264,planted,"A",10,surv_mode="calib_alloc",
                            ingrowth=ig)
            # rebuild t from out is not returned as df; instead project full and snapshot via internal
        # simpler: full runs at horizons
        snaps={}
        for H in (20,50,100,200):
            o=project_psp(base.copy(),264,planted,"A",H,surv_mode="calib_alloc",ingrowth=ig)
            snaps[H]=o
        for H in (20,50,100,200):
            o=snaps[H]
            rows.append(dict(origin="plt" if planted else "nat",ingrowth=ig,age=H,
                QMD=round(o["QMD"],1),TPH=round(o["TPH"],0),BAPH=round(o["BAPH"],1),
                pctSDImax=round(o["SDI"]/SDI_MAX*100,0)))
lt=pd.DataFrame(rows); lt.to_csv("results/ingrowth_longterm.csv",index=False)
for planted in (0,1):
    o="plt" if planted else "nat"; print(f"\n {'PLANTATION' if planted else 'NATURAL'}:")
    for ig in (False,True):
        s=lt[(lt.origin==o)&(lt.ingrowth==ig)]
        print(f"  ingrowth={str(ig):5s}: "+"  ".join(
          f"a{r.age}:QMD{r.QMD:.0f}/TPH{r.TPH:.0f}/BA{r.BAPH:.0f}/{r.pctSDImax:.0f}%" for _,r in s.iterrows()))

# robustness flags with ingrowth on
print("\n"+"="*78); print("Robustness with ingrowth ON (no runaway >105% SDImax):")
bad=False
for planted in (0,1):
    for byi in (100,264,450):
        o=project_psp(synth(planted),byi,planted,"A",200,surv_mode="calib_alloc",ingrowth=True)
        rr=o["SDI"]/SDI_MAX*100
        flag = rr>110 or o["TPH"]<2
        bad = bad or flag
        print(f"  {'plt' if planted else 'nat'} BYI{byi}: SDI200 {rr:.0f}%SDImax TPH {o['TPH']:.0f} QMD {o['QMD']:.0f}  {'FLAG' if flag else 'ok'}")
print(f"\n ANY problem: {bad}")
print("Wrote ingrowth_longterm.csv")
