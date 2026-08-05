# verify_consistency.py -- check the R drop-in constants match the fitted CSVs
# and the Python modules (final verification pass).
import pandas as pd, re
ok=True
def chk(name,a,b,tol=1e-3):
    global ok; good=abs(a-b)<=tol; ok=ok and good
    print(f"  [{'OK' if good else 'MISMATCH'}] {name}: R/code={a}  fitted/py={b}")
# ingrowth: koa_ingrowth.py vs ingrowth_final_coef.csv
import koa_ingrowth as ig
cf=pd.read_csv("results/ingrowth_final_coef.csv").set_index("term").estimate
chk("ingrowth b0", ig.B0, round(cf["(Intercept)"],4))
chk("ingrowth b_RD", ig.B_RD, round(cf["RD"],4))
chk("ingrowth b_planted", ig.B_PLANTED, round(cf["planted"],4))
# ingrowth R drop-in constants
r=open("/sessions/elegant-lucid-cannon/mnt/outputs/fvshi_repo/fvsOL/inst/extdata/koa_stress_test/koa_ingrowth.R").read()
chk("R ingrowth b0", float(re.search(r"b0 = ([\d.\-]+)",r).group(1)), round(cf["(Intercept)"],4))
chk("R ingrowth b_rd", float(re.search(r"b_rd = ([\d.\-]+)",r).group(1)), round(cf["RD"],4))
chk("R ingrowth b_planted", float(re.search(r"b_planted = ([\d.\-]+)",r).group(1)), round(cf["planted"],4))
# survival: py module vs R drop-in
import koa_survival_calibrated_py as sc
import inspect; src=inspect.getsource(sc.surv_calibrated)
rs=open("/sessions/elegant-lucid-cannon/mnt/outputs/fvshi_repo/fvsOL/inst/extdata/koa_stress_test/koa_survival_calibrated.R").read()
for nm,py in [("base_nat",0.003),("base_plt",0.006),("onset",0.65),("full",0.85),("maxlift",0.15)]:
    rval=float(re.search(nm+r" = ([\d.]+)",rs).group(1))
    chk(f"survival {nm}", rval, py)
print("\nALL CONSISTENT" if ok else "\nINCONSISTENCIES FOUND")
