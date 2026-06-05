"""koa_survival_calibrated_py.py
Tuned calibrated koa annual survival. Design follows the data diagnostics
(tune_survival.R): mortality is governed by (1) a low density-independent
background that differs by origin (plantations lower), and (2) a
density-dependent self-thinning term that ramps as relative density (SDI/SDImax)
passes the self-thinning onset. BYI is intentionally NOT a direct mortality
driver: its apparent effect in the data is confined to one cluster of high-BYI
natural plots and cannot be estimated for plantations; BYI influences long-term
density correctly through GROWTH driving stands into self-thinning.
"""
import numpy as np
SDI_MAX = 500.0

def surv_calibrated(dbh, ht, cr, rht, byi, bal=0.0, baph=0.0, planted=0, sdi=None,
                    base_nat=0.005, base_plt=0.003, onset=0.55, ramp=0.55,
                    power=2.0, mort_max=0.15):
    dbh = np.asarray(dbh, float)
    RD = (sdi/SDI_MAX) if sdi is not None else (np.asarray(baph,float)/60.0)
    base = base_plt if planted == 1 else base_nat
    thin = ramp * np.maximum(0.0, RD - onset)**power
    mort = np.clip(base + thin, 0.0, mort_max)
    return 1.0 - mort

def make(**kw):
    def f(dbh, ht, cr, rht, byi, bal=0.0, baph=0.0, planted=0, sdi=None, **k):
        return surv_calibrated(dbh, ht, cr, rht, byi, bal=bal, baph=baph,
                               planted=planted, sdi=sdi, **kw)
    return f
