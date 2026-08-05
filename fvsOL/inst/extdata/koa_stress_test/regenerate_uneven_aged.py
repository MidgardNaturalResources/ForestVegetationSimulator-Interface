"""
regenerate_uneven_aged.py
=========================
Regenerate the uneven-aged natural (ingrowth) scenario for Table 8, Fig 8, and
Section 3.5, from the deposited engine. The original integrated harness was not
preserved; this re-derives the scenario reproducibly from:
  - the component equations (koa_equations, lineage A == HiGy.R),
  - the calibrated density-dependent survival (koa_survival_calibrated_py),
  - the ingrowth model Eq. 8 (koa_ingrowth: exp(5.3836 - 3.0933 RD - 1.6359 planted)),
  - the tree-list engine project_psp with ingrowth on.

Initialization (documented): the uneven-aged natural scenario starts from the same
young natural cohort as the even-aged natural scenario (single class, DBH 8 cm,
500 trees/ha) with annual ingrowth (Eq. 8) enabled. Self-thinning, size caps, and
the HT ceiling are the deployed operational constraints.

Uncertainty: Normal(estimate, SE) draws on height a0/a1/b/c and the ln(BYI) increment
terms (Table S6), plus the survival self-thinning calibration (maxlift). 300 reps.

Outputs: uneven_aged_table8.csv (points + 95% CI), uneven_aged_traj.csv (annual
trajectories for Fig 8), and prints the Section 3.5 summary quantities.
"""
import numpy as np, pandas as pd
import koa_equations as KE
from koa_projector import project_psp, sdi_of, SDI_MAX
from koa_equations import predict_HT
from koa_ingrowth import ingrowth_annual

N_REPS, SEED = 300, 42
AGES = [20, 40, 60, 100]
BYIS = [(100, "Low"), (264, "Medium"), (450, "High")]
HT0 = dict(KE.HT_P); DDBH8 = KE.LineageA.DDBH["b8"]; DHT8 = KE.LineageA.DHT["b8"]
SE = dict(a0=0.614, a1=0.018, b=0.002, c=0.019, dd=0.089, dh=0.081)

def init_stand():
    d0 = 8.0
    return pd.DataFrame(dict(dbh=[d0], ht=[float(predict_HT(d0, 1.0, d0, 264))],
                             cr=[0.6], expf=[500.0]))

def setp(rng):
    if rng is None:
        KE.HT_P.update(HT0); KE.LineageA.DDBH["b8"] = DDBH8; KE.LineageA.DHT["b8"] = DHT8
        return
    KE.HT_P["a0"] = rng.normal(HT0["a0"], SE["a0"]); KE.HT_P["a1"] = rng.normal(HT0["a1"], SE["a1"])
    KE.HT_P["b"] = max(1e-4, rng.normal(HT0["b"], SE["b"])); KE.HT_P["c"] = rng.normal(HT0["c"], SE["c"])
    KE.LineageA.DDBH["b8"] = rng.normal(DDBH8, SE["dd"]); KE.LineageA.DHT["b8"] = rng.normal(DHT8, SE["dh"])

def project(byi, ages):
    """Return dict age-> (QMD, TPH, BAPH, VOL, SDIpct) using return_traj for one run."""
    o = project_psp(init_stand(), byi, 0, "A", max(ages), surv_mode="calib_alloc",
                    ingrowth=True, return_traj=True)
    tj = o["traj"]  # columns BAPH,TPH,QMD,SDI,ingrowth ; index = year-1
    out = {}
    for a in ages:
        if a-1 < len(tj):
            r = tj.iloc[a-1]
            # mean HT not in traj; approximate volume from final-style: use BAPH * Hbar.
            out[a] = dict(QMD=r.QMD, TPH=r.TPH, BAPH=r.BAPH, SDIpct=100*r.SDI/SDI_MAX,
                          ingrowth=r.ingrowth)
        else:
            out[a] = dict(QMD=np.nan, TPH=np.nan, BAPH=np.nan, SDIpct=np.nan, ingrowth=np.nan)
    return out, tj

def hbar_vol(byi, qmd, baph):
    # basal-area-weighted mean height of the average tree, V = BAPH * Hbar * 0.40
    h = float(predict_HT(qmd, baph, qmd, byi))
    return baph * h * 0.40

def main():
    rng = np.random.default_rng(SEED)
    rows = []; trajs = []
    for byi, site in BYIS:
        setp(None); pt, tj = project(byi, AGES)
        tj2 = tj.copy(); tj2["age"] = np.arange(1, len(tj2)+1); tj2["BYI"] = byi; tj2["Site"] = site
        tj2["Vol"] = [hbar_vol(byi, q, b) for q, b in zip(tj2.QMD, tj2.BAPH)]
        tj2["MAI"] = tj2["Vol"]/tj2["age"]
        trajs.append(tj2)
        draws = {a: {"QMD": [], "VOL": []} for a in AGES}
        for _ in range(N_REPS):
            setp(rng); o, _ = project(byi, AGES)
            for a in AGES:
                draws[a]["QMD"].append(o[a]["QMD"])
                draws[a]["VOL"].append(hbar_vol(byi, o[a]["QMD"], o[a]["BAPH"]))
        setp(None)
        for a in AGES:
            vol = hbar_vol(byi, pt[a]["QMD"], pt[a]["BAPH"])
            q = np.array(draws[a]["QMD"]); v = np.array(draws[a]["VOL"])
            rows.append(dict(Scenario="Uneven-aged natural", Site=f"{site} ({byi})", Age=a,
                QMD=round(pt[a]["QMD"],1), QMD_lo=round(np.nanpercentile(q,2.5),1), QMD_hi=round(np.nanpercentile(q,97.5),1),
                TPH=round(pt[a]["TPH"]), BAPH=round(pt[a]["BAPH"],1),
                VOL=round(vol,1), VOL_lo=round(np.nanpercentile(v,2.5),1), VOL_hi=round(np.nanpercentile(v,97.5),1),
                MAI=round(vol/a,2), SDIpct=round(pt[a]["SDIpct"])))
    out = pd.DataFrame(rows); out.to_csv("uneven_aged_table8.csv", index=False)
    pd.concat(trajs).to_csv("uneven_aged_traj.csv", index=False)
    print(out.to_string(index=False))
    # Section 3.5 quantities
    print("\n--- Section 3.5 summary ---")
    for byi, site in BYIS:
        s = out[out.Site.str.startswith(site)]
        tj = [t for t in trajs if t.Site.iloc[0]==site][0]
        qmd100 = s[s.Age==100].QMD.iloc[0]; ba100 = s[s.Age==100].BAPH.iloc[0]
        sdimax = s.SDIpct.max(); ig_young = tj.ingrowth.iloc[:10].mean(); ig_old = tj.ingrowth.iloc[-20:].mean()
        print(f"  {site}: QMD@100={qmd100} cm, BA@100={ba100} m2/ha, peak %SDImax={sdimax}, "
              f"ingrowth young~{ig_young:.0f}, mature~{ig_old:.0f} tph/yr")

if __name__ == "__main__":
    main()
