# koa_survival_calibrated.R
# Recommended operational survival for FVS-HI koa, to REPLACE the published
# cloglog (surv.parm / surv_prob) in HiGy.R.
#
# Rationale: the published per-tree cloglog discriminates well (AUC 0.95) but is
# numerically unstable applied per tree (annual survival 0.80 at CR 0.5, ~0 at
# CR 0.7) and collapses real stands, especially plantations. Free refits on the
# 280-death dataset overfit (separation; flipping origin flips survival). A
# calibrated rate anchored to OBSERVED annual mortality by origin (natural 1.5%,
# plantation 0.5%), modulated by relative density and tree size, is stable,
# monotonic, origin-correct, and tracks observed cohort retention (natural ~53%,
# plantation ~66% at 40 yr). Tune the four constants as more remeasurement data
# accrue.
#
# Inputs (metric): dbh cm, baph m2/ha (plot), planted 0/1.

koa.SURV.calibrated <- function(dbh, baph, planted = 0,
                                base_nat = 0.015, base_plt = 0.005,
                                BA_envelope = 60, thin_onset = 0.55,
                                comp_lift = 2.5, size_lift = 0.8, size_scale = 8,
                                mort_min = 0.002, mort_max = 0.12) {
  base <- ifelse(planted == 1, base_plt, base_nat)
  RD   <- baph / BA_envelope                                  # relative density
  comp <- 1 + comp_lift * pmax(0, RD - thin_onset) / (1 - thin_onset)
  size <- 1 + size_lift * exp(-dbh / size_scale)             # small-tree lift
  mort <- pmin(pmax(base * comp * size, mort_min), mort_max) # annual mortality
  1 - mort                                                    # annual survival
}

# Drop-in for HiGy.R calc_mortality(): replace the surv_prob(...) call with
#   surv = koa.SURV.calibrated(dbh, ba.plot, planted = stand$planted)
# and keep dexpf = expf * (1 - surv) * mort.mult as is.
