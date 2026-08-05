"""koa_ingrowth.py -- koa-only annualized ingrowth as a function of relative
density (RD = SDI/SDImax) and origin. Reconstructed from AK.HT.csv remeasurements
(363 plot-periods) and fit with a quasi-Poisson log model:

    E[ingrowth, trees/ha/yr] = exp( 5.3836 - 3.0933*RD - 1.6359*planted )

RD is the dominant driver (p<1e-4); plantations have ~5x less ingrowth (rate
ratio 0.195, p=0.029). BYI is NOT included: neither a BYI main effect (p=0.28)
nor a BYI x RD interaction (p=0.56) is significant, and % koa BA has no variation
(koa PSPs are pure koa). BYI acts on ingrowth indirectly through growth (it
raises RD faster). An OPTIONAL BYI multiplier (byi_c, default 0 = off) is kept
for when richer data allow testing Aaron's hypothesis that BYI raises ingrowth.
Annualized by construction (rate per year)."""
import numpy as np
SDI_MAX = 500.0
B0, B_RD, B_PLANTED = 5.3836, -3.0933, -1.6359

def ingrowth_annual(sdi=None, baph=None, planted=0, byi=264.0, byi_c=0.0,
                    byi_ref=390.0, rd=None, cap=160.0):
    """Expected annual koa ingrowth (trees/ha/yr) entering at threshold DBH."""
    if rd is None:
        rd = (sdi/SDI_MAX) if sdi is not None else (baph/60.0)
    e = np.exp(B0 + B_RD*rd + B_PLANTED*int(planted))
    if byi_c:
        e = e * (max(byi, 1.0)/byi_ref)**byi_c
    return float(np.clip(e, 0.0, cap))
