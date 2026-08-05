import numpy as np, pandas as pd
from koa_equations import LineageA as A
# Survival sensitivity to CR and rHT at a typical small koa tree
print("Annual survival (LineageA == HiGy.R == FINAL), DBH=12, HT=8, BYI=264")
for cr in [0.3,0.4,0.5,0.6,0.7,0.8]:
    row=[f"CR={cr}"]
    for rht in [0.3,0.5,0.8]:
        s=float(A.surv_annual(12,8,cr,rht,264))
        row.append(f"rHT={rht}:{s:.3f}")
    print("  "+"  ".join(row))
print("\n30-yr cumulative survival at CR,rHT=0.5:")
for cr in [0.4,0.5,0.6,0.7]:
    s=float(A.surv_annual(12,8,cr,0.5,264)); print(f"  CR={cr}: annual {s:.3f} -> 30yr {s**30:.4f}")
# BYI sensitivity
print("\nBYI sensitivity, DBH=12,HT=8,CR=0.5,rHT=0.5:")
for byi in [50,100,200,264,400,600]:
    s=float(A.surv_annual(12,8,0.5,0.5,byi)); print(f"  BYI={byi}: {s:.3f}/yr -> 30yr {s**30:.4f}")

# PSP remeasurement extent
psp=pd.read_csv("/sessions/elegant-lucid-cannon/mnt/Koa/AK_PSP_data_compiled.csv")
print("\n=== PSP extent ===")
print("tree rows:",len(psp),"distinct PSP:",psp.PSP.nunique())
g=psp[psp.DBH_cm>0].groupby("PSP").MeasSeq.nunique()
print("PSP with >=2 measurements:",(g>=2).sum())
# year span per PSP via Yfloor
ys=psp.groupby("PSP").Yfloor.agg(['min','max'])
ys['span']=ys['max']-ys['min']
print("PSP with span>=5 yr:",(ys.span>=5).sum(),"| span>=10:",(ys.span>=10).sum())
print("span distribution:",ys.span.describe()[['min','25%','50%','75%','max']].round(1).to_dict())
