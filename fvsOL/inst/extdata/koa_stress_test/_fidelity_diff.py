import numpy as np, pandas as pd
from koa_equations import predict_HT, predict_HCB, LineageA as A
g = pd.read_csv("results/fidelity_R_HiGy.csv")
py = {}
py["HT"]  = predict_HT(g.dbh.values, g.baph.values, g.qmd.values, g.byi.values)
py["HCB"] = predict_HCB(g.dbh.values, g.ht.values, g.bal.values, g.baph.values, g.byi.values)
py["DDBH"]= A.dDBH(g.dbh.values, g.baph.values, g.bal.values, g.cr.values, g.byi.values, g.planted.values)
py["DHT"] = A.dHT(g.ht.values, g.baph.values, g.bal.values, g.cr.values, g.byi.values, g.planted.values)
py["SURV"]= A.surv_annual(g.dbh.values, g.ht.values, g.cr.values, g.rht.values, g.byi.values)
print("Port fidelity: max|R(HiGy.R) - Python(LineageA)| over %d grid points" % len(g))
for k in ["HT","HCB","DDBH","DHT","SURV"]:
    d = np.abs(py[k]-g[k].values)
    print(f"  {k:5s}  maxdiff={d.max():.3e}  meandiff={d.mean():.3e}")
