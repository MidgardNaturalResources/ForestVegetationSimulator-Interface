# koa_ingrowth.R   (2026-06-05; SDImax corrected 2026-08-05)
# Koa-only annualized ingrowth for FVS-HI, to add the missing recruitment
# component (HiGy.R currently has no ingrowth). Follows the two-stage / density
# logic of Li, Weiskittel & Kershaw (2011, Can. J. For. Res. 41:2077-2089) but
# reduced to a single annualized expectation keyed on relative density (RD), per
# Aaron's preference.
#
# ============================================================================
# 2026-08-05 CHANGE NOTE: SDImax corrected from 500 to 1350, and b_rd rescaled
# so that predictions are unchanged.
#
# WHY SDImax CHANGED. The value 500 was an FIA-only quantity. A Reineke boundary
# fit on 269 pooled plot-years (Cardinal job 13316017), 0.99 quantile with the
# slope fixed at -1.605, gives SDImax = 1352 with a 95% cluster bootstrap
# interval of 1173 to 1882; a stochastic frontier on the same data gives 1087
# (781 to 1344); an FIA-only refit reproduces 708 and identifies the provenance
# of 500. Under SDImax = 500, 41.4% of tree records and 27.9% of plot-years sit
# above RD 1.0 with a maximum RD of 4.39. Under 1350 that falls to 0.95% and
# 0.74% with a maximum RD of 1.63. The deployed value is 1350.
#
# WHY b_rd MOVED AND THE PREDICTIONS DID NOT. Unlike the survival ramp, this is
# a smooth log-linear function of density with no thresholds, so re-expressing
# it under a new SDImax is an EXACT reparameterisation rather than a
# re-tuning. The fitted quantity is a slope per unit of SDI,
#     b_sdi = -3.0933 / 500 = -0.0061866 per SDI unit,
# and the RD-scale coefficient is b_rd = b_sdi * SDImax, which is
# -3.0933 at SDImax 500 and -8.3519 at SDImax 1350. Refitting the quasi-Poisson
# model with RD recomputed at SDImax = 1350 returns exactly this value, because
# the covariate is a linear rescaling of the old one. Verified numerically over
# SDI 0 to 2000 for both origins: maximum absolute difference between the old
# and rescaled predictions is 5.0e-14 trees ha⁻¹ yr⁻¹.
#
# WHAT A NAIVE SWAP WOULD HAVE DONE. Changing SDImax to 1350 while leaving
# b_rd at -3.0933 does not preserve anything. It flattens the density response
# by a factor of 2.7 and inflates recruitment badly: at SDI 500 a natural stand
# would recruit 69.3 instead of 9.9 trees ha⁻¹ yr⁻¹, seven times too many, and
# the 160 trees ha⁻¹ yr⁻¹ cap would bind out to SDI 100 instead of SDI 0. Do
# not do it.
#
# CONSISTENCY REQUIREMENT. koa_survival_calibrated.R and this file must carry
# the SAME SDImax. If one is changed the other must be changed in the same
# commit, or the mortality ramp and the recruitment response will be reading
# two different density scales.
# ============================================================================
#
# DATA. Reconstructed from AK.HT.csv remeasurements (363 plot-periods, 22%
# with ingrowth). Quasi-Poisson log model:
#
#   E[ingrowth, trees ha⁻¹ yr⁻¹] = exp( 5.3836 - 8.3519*RD - 1.6359*planted )
#
#   RD = SDI / SDImax (SDImax = 1350); planted = 0 natural, 1 plantation.
#   Equivalently and invariantly: exp( 5.3836 - 0.0061866*SDI - 1.6359*planted ).
#
# FINDINGS. RD is the dominant driver (p < 1e-4): ingrowth falls from ~160
# (open) to ~13 trees ha⁻¹ yr⁻¹ near canopy closure, which under the corrected
# scale is SDI ~ 340 rather than RD 0.68 of a 500 ceiling. Plantations have ~5x
# less ingrowth (rate ratio 0.195, p = 0.029) -- managed/weeded. BYI is NOT
# included: neither a BYI main effect (p = 0.28) nor a BYI x RD interaction
# (p = 0.56) is significant, and % koa BA has no variation (koa stands are ~pure
# koa). BYI acts on ingrowth indirectly through growth (it raises RD faster).
# byi_c > 0 enables an OPTIONAL, untested BYI multiplier (Aaron's hypothesis:
# BYI raises ingrowth).
#
# Recruits enter at a threshold DBH (default 2.5 cm); set HT from koa.HT and an
# initial crown ratio, then add to the tree list before the next cycle.
#
# NOTE ON PARAMETERISATION. b_rd is now DERIVED from b_sdi and SDImax rather
# than hard-coded, so that any future change to SDImax cannot silently change
# the fitted density response. Passing b_rd directly still works and overrides
# the derivation. The BA fallback reference moves from 60 to 72 m² ha⁻¹ for the
# same reason it does in koa_survival_calibrated.R: 60 m² ha⁻¹ implies an SDImax
# of 927 to 1337 over QMD 20 to 50 cm, which was already incompatible with 500,
# and 72 m² ha⁻¹ is the Reineke-consistent basal area at SDI 1350 and QMD 30 cm.

koa.ingrowth <- function(sdi = NA, baph = NA, planted = 0,
                         SDImax = 1350,             # was 500; Reineke 0.99
                                                    # quantile 1352 (1173-1882)
                         b0 = 5.3836,
                         b_sdi = -0.0061866,        # fitted slope per SDI unit
                                                    # (= -3.0933 / 500)
                         b_rd = b_sdi * SDImax,     # -8.3519 at SDImax 1350
                                                    # (was -3.0933 at SDImax 500)
                         b_planted = -1.6359,
                         byi = 264, byi_c = 0, byi_ref = 390, cap = 160,
                         baph_ref = 72) {           # was 60; see note above
  RD <- if (!is.na(sdi)) sdi / SDImax else baph / baph_ref
  e  <- exp(b0 + b_rd * RD + b_planted * (planted == 1))
  if (byi_c != 0) e <- e * (pmax(byi, 1) / byi_ref)^byi_c
  pmin(pmax(e, 0), cap)                       # expected recruits, trees ha⁻¹ yr⁻¹
}

# Drop-in for HiGy.R (add an ingrowth step, e.g. at FVS stop point 6):
#   n.rec <- koa.ingrowth(sdi = sdi, planted = stand$planted)   # trees ha⁻¹ yr⁻¹
#   if (n.rec > 0.01) add a recruit record: dbh = 2.5 cm,
#       ht = koa.HT(2.5, baph, qmd, byi), cr ~ 0.6, expf = n.rec
# Verified on 200-yr projections: stands become realistically multi-cohort
# (sustained BA/density rather than self-thinning to a few large trees), with no
# runaway or collapse across origin and BYI. Re-verified 2026-08-05 to 300 years
# against the re-tuned survival ramp at SDImax 1350: with this ingrowth active
# the projection reaches a genuine steady state, natural 16.3 m² ha⁻¹ at SDI 367
# and plantation 15.9 m² ha⁻¹ at SDI 301, with density drift under 0.02% per 20
# years. Without ingrowth the same runs are still drifting at 6 to 13% per 20
# years at year 300 and settle 35% lower in basal area, so this component is
# what actually closes the long-run stand dynamics; wiring it into HiGy.R is a
# prerequisite for enabling the recovered mortality levels documented in
# koa_survival_calibrated.R.
