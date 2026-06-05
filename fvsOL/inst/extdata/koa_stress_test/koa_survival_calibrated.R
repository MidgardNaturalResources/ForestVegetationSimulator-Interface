# koa_survival_calibrated.R   (tuned 2026-06-05)
# Recommended operational survival for FVS-HI koa, to REPLACE the published
# cloglog (surv.parm / surv_prob) in HiGy.R.
#
# WHY NOT A FITTED GLM. The published per-tree cloglog discriminates well
# (AUC 0.95) but is numerically unstable applied per tree (annual survival 0.80
# at CR 0.5, ~0 at CR 0.7) and collapses real stands, especially plantations.
# Data diagnostics (tune_survival.R) show the survival signal cannot support a
# free per-tree GLM: only 280 deaths; mortality is lowest in small trees
# (0-5 cm: 0.04%/yr) not highest; and the apparent BYI effect is an artifact of
# ONE cluster of high-BYI natural plots (BYI>408: 3.3%/yr vs ~0.15% otherwise),
# a cluster that contains no plantations, so a BYI mortality effect cannot even
# be estimated for plantations. Forcing BYI into a GLM yields absurd, unstable
# rate ratios (~890x per log unit).
#
# DESIGN. Annual mortality = a low density-independent background that differs
# by origin (plantations lower, as observed) PLUS a density-dependent
# self-thinning term that starts at 0.65 relative density (SDI/SDImax), increases
# LINEARLY to a maximum lift at 0.85 RD, and plateaus above. BYI is deliberately
# NOT a direct mortality driver; it influences long-term density correctly
# through GROWTH (higher BYI reaches the self-thinning boundary sooner). Verified
# on 200-year projections across origin and BYI: SDI peaks at 72 to 88% of
# SDImax then tracks the self-thinning line (no runaway, no collapse, QMD
# monotonic); natural QMD approaches ~90 cm, plantations plateau at the 60 cm cap.
#
# Inputs (metric): sdi (stand SDI), baph m2/ha (fallback if sdi missing),
# planted 0/1.

koa.SURV.calibrated <- function(sdi = NA, baph = NA, planted = 0,
                                base_nat = 0.005, base_plt = 0.003,
                                SDImax = 500, onset = 0.65, full = 0.85,
                                maxlift = 0.15, mort_max = 0.20) {
  RD   <- if (!is.na(sdi)) sdi / SDImax else baph / 60          # relative density
  base <- ifelse(planted == 1, base_plt, base_nat)
  frac <- pmin(pmax((RD - onset) / (full - onset), 0), 1)       # 0 at .65, 1 at .85+
  mort <- pmin(pmax(base + maxlift * frac, 0), mort_max)        # annual mortality
  1 - mort                                                      # annual survival
}

# Drop-in for HiGy.R calc_mortality(): compute stand SDI from the plot summary
#   (sdi = tph.plot * (qmd/25)^1.6), then
#   surv = koa.SURV.calibrated(sdi = sdi, planted = stand$planted)
# and keep dexpf = expf * (1 - surv) * mort.mult as is. The mortality is a stand
# (density) rate applied to all trees; refine the four constants as Kahikinui,
# KMR, and Kualoa remeasurements accrue.
