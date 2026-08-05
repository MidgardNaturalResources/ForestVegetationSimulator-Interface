"""
regenerate_fig8.py
==================
Figure 8: koa stand projections for three management scenarios (columns:
even-aged natural, even-aged planted, uneven-aged natural with ingrowth) across
Low/Medium/High BYI (colors), three metrics (rows: QMD, stem volume, MAI),
with 95% Monte Carlo confidence ribbons.

Even-aged scenarios: project_cohort (the engine behind Table 8 even-aged).
Uneven-aged: project_psp with ingrowth (Eq. 8), initialized as the even-aged
natural cohort. Perturbation: Normal(estimate, SE) on height a0/a1/b/c and the
ln(BYI) increment terms, plus the survival self-thinning calibration. 200 reps.

Run: python regenerate_fig8.py   (writes Fig8.png, 320 dpi)
"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import koa_equations as KE
from koa_projector import project_cohort, project_psp
from koa_equations import predict_HT
from koa_survival_calibrated_py import make

N, SEED, MAXA = 200, 42, 100
BYIS = [(100, "Low"), (264, "Medium"), (450, "High")]
COL = {"Low": "#2166ac", "Medium": "#1a9850", "High": "#d73027"}
HT0 = dict(KE.HT_P); DD8 = KE.LineageA.DDBH["b8"]; DH8 = KE.LineageA.DHT["b8"]
SE = dict(a0=0.614, a1=0.018, b=0.002, c=0.019, dd=0.089, dh=0.081)
ages = np.arange(1, MAXA + 1)

def setp(rng):
    if rng is None:
        KE.HT_P.update(HT0); KE.LineageA.DDBH["b8"] = DD8; KE.LineageA.DHT["b8"] = DH8; return 0.15
    KE.HT_P["a0"] = rng.normal(HT0["a0"], SE["a0"]); KE.HT_P["a1"] = rng.normal(HT0["a1"], SE["a1"])
    KE.HT_P["b"] = max(1e-4, rng.normal(HT0["b"], SE["b"])); KE.HT_P["c"] = rng.normal(HT0["c"], SE["c"])
    KE.LineageA.DDBH["b8"] = rng.normal(DD8, SE["dd"]); KE.LineageA.DHT["b8"] = rng.normal(DH8, SE["dh"])
    return float(np.clip(rng.normal(0.15, 0.03), 0.05, 0.25))

def hbar_vol(byi, qmd, baph):
    return baph * float(predict_HT(qmd, max(baph, 0.1), max(qmd, 1.0), byi)) * 0.40

def run_even(byi, planted, ml):
    df = project_cohort(byi, bool(planted), "A", bounded=True, surv_fn=make(maxlift=ml), max_age=MAXA)
    df = df.set_index("age").reindex(ages)
    qmd = df["QMD"].to_numpy(); vol = (df["BAPH"]*df["HT"]*0.40).to_numpy()
    return qmd, vol, vol/ages

def run_uneven(byi, ml):
    init = pd.DataFrame(dict(dbh=[8.0], ht=[float(predict_HT(8.0,1.0,8.0,264))], cr=[0.6], expf=[500.0]))
    o = project_psp(init, byi, 0, "A", MAXA, surv_mode="calib_alloc", ingrowth=True, return_traj=True)
    tj = o["traj"]
    qmd = np.full(MAXA, np.nan); vol = np.full(MAXA, np.nan)
    for i in range(min(MAXA, len(tj))):
        qmd[i] = tj.iloc[i].QMD; vol[i] = hbar_vol(byi, tj.iloc[i].QMD, tj.iloc[i].BAPH)
    return qmd, vol, vol/ages

SCEN = [("Even-aged natural", lambda byi, ml: run_even(byi, 0, ml)),
        ("Even-aged planted", lambda byi, ml: run_even(byi, 1, ml)),
        ("Uneven-aged natural", lambda byi, ml: run_uneven(byi, ml))]
METRICS = ["QMD (cm)", "Stem volume (m³/ha)", "MAI (m³/ha/yr)"]

def collect():
    rng = np.random.default_rng(SEED)
    data = {}  # (scen, site) -> dict metric -> (mean, lo, hi)
    for sname, fn in SCEN:
        for byi, site in BYIS:
            ml = setp(None); base = fn(byi, ml)
            stacks = [np.empty((N, MAXA)) for _ in range(3)]
            for r in range(N):
                ml = setp(rng); res = fn(byi, ml)
                for m in range(3): stacks[m][r] = res[m]
            setp(None)
            d = {}
            for m in range(3):
                d[m] = (base[m],
                        np.nanpercentile(stacks[m], 2.5, axis=0),
                        np.nanpercentile(stacks[m], 97.5, axis=0))
            data[(sname, site)] = d
    return data

def main():
    data = collect()
    fig, ax = plt.subplots(3, 3, figsize=(12, 9), sharex=True)
    for j, (sname, _) in enumerate(SCEN):
        for i in range(3):
            a = ax[i, j]
            for byi, site in BYIS:
                mean, lo, hi = data[(sname, site)][i]
                a.fill_between(ages, lo, hi, color=COL[site], alpha=0.15, linewidth=0)
                a.plot(ages, mean, color=COL[site], lw=1.6, label=site)
            if i == 0: a.set_title(sname, fontsize=11, fontweight="bold")
            if j == 0: a.set_ylabel(METRICS[i], fontsize=10)
            if i == 2: a.set_xlabel("Stand age (yr)", fontsize=10)
            a.spines[["top", "right"]].set_visible(False)
    ax[0, 0].legend(title="Site (BYI)", frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig("Fig8.png", dpi=320, bbox_inches="tight")
    print("wrote Fig8.png")

if __name__ == "__main__":
    main()
