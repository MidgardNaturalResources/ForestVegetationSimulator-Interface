# koa_ingrowth.R   (2026-06-05)
# Koa-only annualized ingrowth for FVS-HI, to add the missing recruitment
# component (HiGy.R currently has no ingrowth). Follows the two-stage / density
# logic of Li, Weiskittel & Kershaw (2011, Can. J. For. Res. 41:2077-2089) but
# reduced to a single annualized expectation keyed on relative density (RD), per
# Aaron's preference.
#
# DATA. Reconstructed from AK.HT.csv remeasurements (363 plot-periods, 22%
# with ingrowth). Quasi-Poisson log model:
#
#   E[ingrowth, trees/ha/yr] = exp( 5.3836 - 3.0933*RD - 1.6359*planted )
#
#   RD = SDI / SDImax (SDImax = 500); planted = 0 natural, 1 plantation.
#
# FINDINGS. RD is the dominant driver (p < 1e-4): ingrowth falls from ~160
# (open) to ~13 trees/ha/yr near canopy closure. Plantations have ~5x less
# ingrowth (rate ratio 0.195, p = 0.029) -- managed/weeded. BYI is NOT included:
# neither a BYI main effect (p = 0.28) nor a BYI x RD interaction (p = 0.56) is
# significant, and % koa BA has no variation (koa stands are ~pure koa). BYI acts
# on ingrowth indirectly through growth (it raises RD faster). byi_c > 0 enables
# an OPTIONAL, untested BYI multiplier (Aaron's hypothesis: BYI raises ingrowth).
#
# Recruits enter at a threshold DBH (default 2.5 cm); set HT from koa.HT and an
# initial crown ratio, then add to the tree list before the next cycle.

koa.ingrowth <- function(sdi = NA, baph = NA, planted = 0, SDImax = 500,
                         b0 = 5.3836, b_rd = -3.0933, b_planted = -1.6359,
                         byi = 264, byi_c = 0, byi_ref = 390, cap = 160) {
  RD <- if (!is.na(sdi)) sdi / SDImax else baph / 60
  e  <- exp(b0 + b_rd * RD + b_planted * (planted == 1))
  if (byi_c != 0) e <- e * (pmax(byi, 1) / byi_ref)^byi_c
  pmin(pmax(e, 0), cap)                       # expected recruits, trees/ha/yr
}

# Drop-in for HiGy.R (add an ingrowth step, e.g. at FVS stop point 6):
#   n.rec <- koa.ingrowth(sdi = sdi, planted = stand$planted)   # trees/ha/yr
#   if (n.rec > 0.01) add a recruit record: dbh = 2.5 cm,
#       ht = koa.HT(2.5, baph, qmd, byi), cr ~ 0.6, expf = n.rec
# Verified on 200-yr projections: stands become realistically multi-cohort
# (sustained BA/density rather than self-thinning to a few large trees), with no
# runaway or collapse across origin and BYI.
