"""
regenerate_table8_CI.py
=======================
Regenerate the Table 8 stand-projection means and 95% confidence intervals from
the CALIBRATED deployed engine, so the point estimates and the uncertainty bands
come from the same model.

Engine (point estimates): project_cohort(byi, planted, "A", bounded=True,
        surv_fn=make(maxlift=0.15)) from koa_projector + koa_equations +
        koa_survival_calibrated_py. This reproduces the v45 Table 8 means.

Uncertainty: parameters perturbed with Normal(estimate, SE) draws (Table S6):
  height a0, a1, b, c; the ln(BYI) terms of the diameter and height increment
  models (the BYI sensitivity that separates the site classes). The projection
  survival is the calibrated self-thinning form (not the fitted Table 6 cloglog),
  so its uncertainty is propagated through the self-thinning lift parameter
  (maxlift ~ Normal(0.15, 0.03), clipped), which governs density and therefore
  stocking and volume.

Run: python regenerate_table8_CI.py   (writes Table8_CI_recomputed.csv)
"""
import numpy as np, pandas as pd
import koa_equations as KE
from koa_projector import project_cohort
from koa_survival_calibrated_py import make

N_REPS = 500
SEED = 42
AGES = [20, 40, 60, 100]
SCENARIOS = [("Even-aged natural", 0, 100, "Low"),
             ("Even-aged natural", 0, 264, "Medium"),
             ("Even-aged natural", 0, 450, "High"),
             ("Even-aged planted", 1, 100, "Low"),
             ("Even-aged planted", 1, 264, "Medium"),
             ("Even-aged planted", 1, 450, "High")]

# Table S6 standard errors (same source the manuscript MC used)
SE = dict(ht_a0=0.614, ht_a1=0.018, ht_b=0.002, ht_c=0.019,
          ddbh_lnBYI=0.089, dht_lnBYI=0.081)
MAXLIFT, MAXLIFT_SE = 0.15, 0.03

# baseline (unperturbed) parameter values
HT0 = dict(KE.HT_P)
DDBH_b8_0 = KE.LineageA.DDBH["b8"]   # ln(BYI) term, diameter increment
DHT_b8_0  = KE.LineageA.DHT["b8"]    # ln(BYI) term, height increment

VARS = ["QMD", "HT", "BAPH", "TPH", "VOL", "MAI"]

def run_one(byi, planted, maxlift):
    df = project_cohort(byi, bool(planted), "A", bounded=True,
                        surv_fn=make(maxlift=maxlift), max_age=100)
    df = df.set_index("age")
    out = {}
    for a in AGES:
        if a in df.index:
            r = df.loc[a]
            vol = r["BAPH"] * r["HT"] * 0.40
            out[a] = dict(QMD=r["QMD"], HT=r["HT"], BAPH=r["BAPH"], TPH=r["TPH"],
                          VOL=vol, MAI=vol / a)
        else:
            out[a] = {v: np.nan for v in VARS}
    return out

def set_params(rng=None):
    """Set KE parameters to baseline (rng=None) or a perturbed draw."""
    if rng is None:
        KE.HT_P.update(HT0)
        KE.LineageA.DDBH["b8"] = DDBH_b8_0
        KE.LineageA.DHT["b8"] = DHT_b8_0
        return MAXLIFT
    KE.HT_P["a0"] = rng.normal(HT0["a0"], SE["ht_a0"])
    KE.HT_P["a1"] = rng.normal(HT0["a1"], SE["ht_a1"])
    KE.HT_P["b"]  = max(1e-4, rng.normal(HT0["b"],  SE["ht_b"]))
    KE.HT_P["c"]  = rng.normal(HT0["c"],  SE["ht_c"])
    KE.LineageA.DDBH["b8"] = rng.normal(DDBH_b8_0, SE["ddbh_lnBYI"])
    KE.LineageA.DHT["b8"]  = rng.normal(DHT_b8_0,  SE["dht_lnBYI"])
    return float(np.clip(rng.normal(MAXLIFT, MAXLIFT_SE), 0.05, 0.25))

def main():
    rng = np.random.default_rng(SEED)
    rows = []
    for scen, planted, byi, site in SCENARIOS:
        # point estimate (baseline parameters)
        ml = set_params(None)
        pt = run_one(byi, planted, ml)
        # Monte Carlo
        draws = {a: {v: [] for v in VARS} for a in AGES}
        for _ in range(N_REPS):
            ml = set_params(rng)
            o = run_one(byi, planted, ml)
            for a in AGES:
                for v in VARS:
                    draws[a][v].append(o[a][v])
        set_params(None)  # restore baseline
        for a in AGES:
            row = dict(Scenario=scen, Site=f"{site} ({byi})", Age=a)
            for v in VARS:
                arr = np.array(draws[a][v], float)
                row[v] = round(float(pt[a][v]), 2)
                row[f"{v}_lo"] = round(float(np.nanpercentile(arr, 2.5)), 2)
                row[f"{v}_hi"] = round(float(np.nanpercentile(arr, 97.5)), 2)
            rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv("Table8_CI_recomputed.csv", index=False)
    # readable summary
    pd.set_option("display.width", 200)
    for v in ["QMD", "VOL", "MAI"]:
        print(f"\n=== {v} (point [95% CI]) ===")
        for _, r in out.iterrows():
            print(f"  {r.Scenario[:16]:<16} {r.Site:<12} age {r.Age:>3}: "
                  f"{r[v]:>7} [{r[v+'_lo']:>7}, {r[v+'_hi']:>7}]")
    print("\nwrote Table8_CI_recomputed.csv")

if __name__ == "__main__":
    main()
