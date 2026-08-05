"""run_alloc_test.py -- verify the individual-tree mortality allocation is
constrained to the stand-level rate, and check its long-term effect (preferential
removal of suppressed trees -> faster QMD, same density trajectory)."""
import numpy as np, pandas as pd
from koa_projector import project_psp, sdi_of, SDI_MAX
from koa_survival_calibrated_py import surv_calibrated
from koa_equations import predict_HT
np.random.seed(3)

def synth(planted, n=40):
    dbh=np.clip(np.random.normal(10 if planted else 14,4,n),2,40)
    tph=1400 if planted else 700
    return pd.DataFrame(dict(dbh=dbh, ht=predict_HT(dbh,15.,np.sqrt((dbh**2).mean()),264),
                             expf=np.full(n,tph/n), cr=0.6))

# 1. CONSTRAINT CHECK: one annual step, does expf-weighted mortality == stand rate?
print("="*74); print("CONSTRAINT CHECK (1 yr): realized stand mortality vs calibrated target")
for planted in (0,1):
    t=synth(planted)
    tph0=t.expf.sum()
    out=project_psp(t.copy(),264,planted,"A",1,surv_mode="calib_alloc")
    realized=1-out["TPH"]/tph0
    ba=(t.dbh**2*0.00007854*t.expf).sum(); qmd=np.sqrt(ba/(0.00007854*tph0)); sdi=sdi_of(tph0,qmd)
    target=1-float(surv_calibrated(qmd,15,.5,.5,264,sdi=sdi,planted=planted))
    print(f"  {'plt' if planted else 'nat':3s} RD={sdi/SDI_MAX:.2f}: target {target:.4f}/yr  realized {realized:.4f}/yr")

# 2. LONG-TERM: allocation vs uniform stand rate (QMD + density)
print("\n"+"="*74); print("LONG-TERM (100 yr): uniform stand rate vs size-allocated")
for planted in (0,1):
    t0=synth(planted)
    # uniform: surv_fn that returns scalar stand survival for all trees
    def uni(dbh,ht,cr,rht,byi,bal=0,baph=0,planted=planted,sdi=None,**k):
        s=float(surv_calibrated(np.sqrt(np.mean(dbh**2)),ht,cr,rht,byi,baph=baph,planted=planted,sdi=sdi))
        return np.full(np.shape(dbh), s)
    u=project_psp(t0.copy(),264,planted,"A",100,surv_fn=uni)
    a=project_psp(t0.copy(),264,planted,"A",100,surv_mode="calib_alloc")
    print(f"  {'PLANTATION' if planted else 'NATURAL':10s}: "
          f"uniform  QMD {u['QMD']:.1f} TPH {u['TPH']:.0f} SDI {u['SDI']:.0f} | "
          f"allocated QMD {a['QMD']:.1f} TPH {a['TPH']:.0f} SDI {a['SDI']:.0f}")
print("\n(allocation should give similar TPH/SDI but higher QMD: suppressed trees removed first)")
