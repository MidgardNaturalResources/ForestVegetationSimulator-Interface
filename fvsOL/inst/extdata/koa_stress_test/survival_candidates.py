"""survival_candidates.py -- generic evaluator for fitted survival variants
read from results/survival_candidates.csv. Returns ANNUAL P(alive)."""
import numpy as np, pandas as pd, os
_CSV = os.path.join(os.path.dirname(__file__), "results", "survival_candidates.csv")

def _load():
    # term field can contain commas (e.g. pmax(BAL.0, 0)); parse: first 2 fields
    # are variant,link; last is estimate; middle (joined) is the term.
    out = {}
    with open(_CSV) as fh:
        next(fh)
        for line in fh:
            line = line.rstrip("\n")
            if not line: continue
            parts = line.split(",")
            variant, link = parts[0], parts[1]
            est = float(parts[-1]); term = ",".join(parts[2:-1])
            out.setdefault(variant, dict(link=link, terms={}))["terms"][term] = est
    return out
VARIANTS = _load()

def _feat(term, d, h, cr, rht, byi, bal, baph, pl):
    L = np.log
    d=np.maximum(d,0.1); h=np.maximum(h,0.1); cr=np.clip(cr,0.01,0.999); byi=np.maximum(byi,1.0)
    bal=np.maximum(bal,0.0); baph=np.maximum(baph,0.0)
    table = {
        "(Intercept)": np.ones_like(d*1.0),
        "HT.0": h, "I(log(HT.0))": L(h), "rHT.0": rht, "I(log(CR.0))": L(cr),
        "I(log(HT.0/(DBH.0/100)))": L(h/(d/100)),
        "I(log(BYI/100))": L(byi/100), "I(BYI/1000)": byi/1000, "log(BYI)": L(byi),
        "I(log(DBH.0))": L(d), "DBH.0": d,
        "I((BAL.0 + 1)/log(DBH.0 + 1))": (bal+1)/L(d+1),
        "Planted": np.full_like(d*1.0, float(pl)) if np.isscalar(pl) else pl,
        "I(log(pmax(BAL.0, 0) + 1))": L(bal+1),
        "I(sqrt(BAPH.0 * DBH.0))": np.sqrt(baph*d),
    }
    return table[term]

def predict(variant, dbh, ht, cr, rht, byi, bal=0.0, baph=0.0, planted=0):
    spec = VARIANTS[variant]
    dbh=np.asarray(dbh,float); ht=np.asarray(ht,float); cr=np.asarray(cr,float)
    rht=np.asarray(rht,float)
    eta = np.zeros_like(dbh*1.0)
    for term, est in spec["terms"].items():
        eta = eta + est * _feat(term, dbh, ht, cr, rht, byi, bal, baph, planted)
    if spec["link"] == "cloglog":
        p = 1 - np.exp(-np.exp(np.clip(eta, -50, 50)))
    else:
        p = 1/(1+np.exp(-np.clip(eta, -50, 50)))
    return np.clip(p, 0, 1)

def make_surv_fn(variant):
    def f(dbh, ht, cr, rht, byi, bal=0.0, baph=0.0, planted=0, **k):
        return predict(variant, dbh, ht, cr, rht, byi, bal=bal, baph=baph, planted=planted)
    return f
