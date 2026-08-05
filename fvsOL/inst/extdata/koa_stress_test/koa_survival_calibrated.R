# koa_survival_calibrated.R   (tuned 2026-06-05; SDImax corrected and ramp
#                              re-tuned 2026-08-05)
# Recommended operational survival for FVS-HI koa, to REPLACE the published
# cloglog (surv.parm / surv_prob) in HiGy.R.
#
# ============================================================================
# 2026-08-05 CHANGE NOTE: SDImax corrected from 500 to 1350, and the
# self-thinning ramp RE-TUNED against observed mortality at the corrected scale.
#
# WHY SDImax CHANGED. The value 500 was an FIA-only quantity carried into a
# pooled-sample equation. A Reineke boundary fit on 269 plot-years pooling all
# four sources (Cardinal job 13316017), 0.99 quantile regression with the slope
# fixed at -1.605, gives SDImax = 1352 with a 95% cluster bootstrap interval of
# 1173 to 1882; a normal half-normal stochastic frontier on the same data gives
# 1087 (781 to 1344); refitting on FIA plot-years alone reproduces 708, which
# identifies the provenance of the published 500. Under SDImax = 500, 41.4% of
# tree records and 27.9% of plot-years sit above relative density 1.0, with a
# maximum RD of 4.39, which is impossible by construction. Under SDImax = 1350
# those fractions fall to 0.95% and 0.74% with a maximum RD of 1.63. The
# deployed point value is 1350.
#
# An independent internal check points the same way. The BA fallback in this
# function treats RD = BA / 60 m² ha⁻¹. Setting that equal to the Reineke path,
# SDImax = 60 / (0.013540 * QMD^0.4), gives an implied SDImax of 927 to 1337
# over QMD 20 to 50 cm. The published equation therefore already contained two
# mutually inconsistent density scales that differed by a factor of about 2.3 on
# a typical koa stand. At SDImax = 1350 with the BA reference moved to 72 m²
# ha⁻¹ the two paths agree to within 2%.
#
# WHY THE RAMP WAS RE-TUNED RATHER THAN RESCALED. The thresholds 0.65 and 0.85
# were themselves fitted quantities, tuned against observed mortality UNDER
# SDImax = 500, which placed the self-thinning onset near SDI 325 and full lift
# near SDI 425. Two options were evaluated and the second was adopted:
#   (a) preserve the absolute thresholds 325 and 425 and re-express them as
#       RD 0.2407 and 0.3148 at SDImax = 1350. Dynamically identical to the
#       published equation; only the reporting denominator changes.
#   (b) re-tune the ramp from the data at the corrected scale. ADOPTED.
# Both were tested by profile likelihood on the recovered variant (ii) survival
# data restricted to DBH >= 2.5 cm, the tree list FVS actually carries (n =
# 4141 records, 927 events, 15432 tree-years, 48 installation clusters), with
# base_nat, base_plt and maxlift maximised out at every threshold pair so that
# threshold location is not confounded with mortality level. Option (a) is
# rejected, likelihood ratio 67.3 on 2 df, P = 2.4e-15. Simply keeping 0.65 and
# 0.85 and swapping SDImax to 1350 is rejected far harder, likelihood ratio
# 137.2, P = 1.6e-30, and in that fit maxlift collapses to 1.7e-11, that is, the
# self-thinning term switches off entirely and the model degenerates to a
# constant hazard. That is the failure mode a naive rescale would have shipped.
#
# THE RE-TUNED THRESHOLDS, in absolute SDI so that they do not depend on the
# disputed constant: onset at SDI 200 (95% profile region 75 to 225) and full
# lift at SDI 850 (95% profile region 800 to 925). At SDImax = 1350 that is
# onset RD 0.1481 (0.0556 to 0.1667) and full RD 0.6296 (0.5926 to 0.6852). The
# full-lift threshold is the well identified one and is stable across samples;
# the onset is not, and moves to SDI 750 if sub-2.5 cm seedling records are left
# in. [UNKNOWN: the onset location is sample dependent and should be revisited
# when the sub-2.5 cm records are adjudicated.]
#
# WHAT WAS DELIBERATELY NOT CHANGED. base_nat, base_plt and maxlift are LEVELS,
# and the level of koa mortality is currently controlled by an unresolved
# death-recording convention (the sentinel rows) rather than by SDImax. Fitting
# them on the recovered data gives base_nat 0.0384 (0.0142 to 0.0510), base_plt
# 0.0142 (0.0000 to 0.0218) and maxlift 0.0870 (0.0393 to 0.1436), all 95%
# cluster bootstrap over installation. Those values also REVERSE the published
# origin ordering: observed annual mortality is 6.41% natural against 2.32%
# planted on this sample, and the independent cloglog respecification finds a
# Planted coefficient of -1.151 (wild cluster bootstrap P < 0.001), so planted
# stock is protective, not the reverse. They are recorded here as
# base_nat_recovered, base_plt_recovered and maxlift_recovered but are NOT the
# defaults, because HiGy.R as it stands carries no ingrowth, and with no
# recruitment a 3.8% per year background draws an operational stand down to 2 to
# 3 m² ha⁻¹ of basal area by year 200. With koa_ingrowth active they are safe
# (steady state 13.5 to 16.3 m² ha⁻¹ at 300 years). Enable them only once the
# sentinel convention is settled AND ingrowth is wired into HiGy.R.
# ============================================================================
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
# by origin PLUS a density-dependent self-thinning term that starts at SDI 200
# (RD 0.148 at SDImax 1350), increases LINEARLY to a maximum lift at SDI 850
# (RD 0.630), and plateaus above. BYI is deliberately NOT a direct mortality
# driver; it influences long-term density correctly through GROWTH (higher BYI
# reaches the self-thinning boundary sooner). Re-verified 2026-08-05 on
# 300-year projections, 5 extreme starting states x 2 origins x 3 BYI levels =
# 30 scenarios, plus the original harness starting states: no NaN, no collapse,
# basal area bounded, monotone diameter in every scenario except a plantation
# started at 80 cm against the 60 cm plantation size cap (an input-validation
# artifact, not an equation failure). With koa_ingrowth active the projection
# reaches a genuine steady state, natural 16.3 m² ha⁻¹ at SDI 367 (27.2% of
# SDImax) and plantation 15.9 m² ha⁻¹ at SDI 301 (22.3%), against 15.5 and 18.3
# m² ha⁻¹ under the published equation. Without ingrowth no configuration,
# published or re-tuned, reaches steady state within 300 years; the published
# equation is still drifting at 6 to 13% per 20 years at year 300, so the
# earlier claim that all stands reach steady state was a 200-year artifact and
# does not hold at 300 years. Peak density is 26 to 44% of the corrected SDImax,
# not the 72 to 88% of the old one, which is the same absolute peak expressed
# against a scale that is 2.7 times larger.
#
# Inputs (metric): sdi (stand SDI), baph m² ha⁻¹ (fallback if sdi missing),
# planted 0/1.

# ORIGIN AND BYI (data-checked). Origin: yes, but the SIGN is now disputed. The
# original diagnostic chain showed plantations HIGHER (0.66%/yr) than natural
# (0.22%/yr), read as young-stand establishment mortality, which is why
# base_plt > base_nat below. The recovered mortality data reverse it (natural
# 6.41%/yr against planted 2.32%/yr). The defaults retain the published ordering
# because they retain the published levels; see the change note above. The
# difference is 0.3 percentage points and is swamped by the ramp in any stand
# above SDI 200. BYI: no direct term. Plantations occur only at BYI <= 191
# (max 191; zero records above 399), so an origin x BYI interaction is not
# identifiable (confounded). The real dynamics you would expect -- plantations
# and higher-BYI sites self-thin faster, and survival is lower at higher BYI --
# EMERGE from growth driving stands into the self-thinning ramp sooner
# (verified: natural self-thinning onset age 14 -> 7 as BYI rises 100 -> 550;
# plantations onset 5-8 yr, earlier at every BYI). No BYI coefficient.
#
# NOTE ON PARAMETERISATION. onset and full are now DERIVED from absolute SDI
# thresholds and SDImax rather than hard-coded on the RD scale. That is the
# whole lesson of the 2026-08-05 correction: the ramp is a property of stand
# density, not of the constant chosen to normalise it, so changing SDImax must
# move the RD thresholds and leave the absolute ones alone. Passing onset or
# full directly still works and overrides the derivation.
koa.SURV.calibrated <- function(sdi = NA, baph = NA, planted = 0,
                                base_nat = 0.003, base_plt = 0.006,
                                SDImax = 1350,          # was 500; Reineke 0.99
                                                        # quantile 1352 (1173-1882)
                                onset_sdi = 200,        # was 325 (= 0.65*500);
                                                        # profile 75-225
                                full_sdi  = 850,        # was 425 (= 0.85*500);
                                                        # profile 800-925
                                onset = onset_sdi/SDImax,   # 0.1481 at SDImax 1350
                                full  = full_sdi/SDImax,    # 0.6296 at SDImax 1350
                                maxlift = 0.15, mort_max = 0.20,
                                baph_ref = 72) {        # was 60; Reineke-consistent
                                                        # BA at SDI 1350, QMD 30 cm
  RD   <- if (!is.na(sdi)) sdi / SDImax else baph / baph_ref    # relative density
  base <- ifelse(planted == 1, base_plt, base_nat)
  frac <- pmin(pmax((RD - onset) / (full - onset), 0), 1)  # 0 at SDI 200, 1 at 850+
  mort <- pmin(pmax(base + maxlift * frac, 0), mort_max)   # annual mortality
  1 - mort                                                 # annual survival
}

# ---- Re-estimated levels, NOT the defaults ----------------------------------
# Maximum likelihood on the recovered variant (ii) survival data, DBH >= 2.5 cm,
# with 95% cluster bootstrap intervals over installation (435 of 500 replicates
# retained). Realised annual mortality on the observed record under these values
# is 8.07% natural and 2.37% planted, against observed 6.41% and 2.32%; under
# the retained defaults it is 7.59% and 2.25%. Do not enable until the sentinel
# death-recording convention is resolved and koa_ingrowth is active in HiGy.R.
koa.SURV.levels_recovered <- list(
  base_nat_recovered = 0.0384,   # 95% CI 0.0142 to 0.0510  (published 0.003)
  base_plt_recovered = 0.0142,   # 95% CI 0.0000 to 0.0218  (published 0.006)
  maxlift_recovered  = 0.0870    # 95% CI 0.0393 to 0.1436  (published 0.15)
)

# ---- Allocate the stand mortality rate to individual trees -------------------
# Distributes the calibrated STAND annual mortality across trees by relative
# size (suppressed trees die first), CONSTRAINED so the expansion-factor-weighted
# mean per-tree mortality equals the stand rate. This keeps the validated
# stand-level density trajectory while making self-thinning size-realistic
# (verified: realized stand mortality matches target exactly; over 100 yr it
# preserves SDI/TPH but raises natural QMD ~6 cm by removing small trees).
# Unaffected by the SDImax correction: it takes the stand rate as given.
#
#   dbh    tree DBH (cm); qmd stand QMD (cm); rDBH = dbh/qmd
#   expf   tree expansion factor (trees ha⁻¹)
#   m_stand = 1 - koa.SURV.calibrated(sdi = sdi, planted = planted)
#   beta   concentration of mortality on small trees (default 3)

koa.SURV.allocate <- function(dbh, expf, qmd, m_stand, beta = 3) {
  rDBH <- dbh / pmax(qmd, 0.1)
  w    <- exp(-beta * (rDBH - 1))                 # small trees (rDBH<1) -> w>1
  wbar <- sum(w * expf) / sum(expf)               # expf-weighted mean weight
  pmin(pmax(m_stand * w / wbar, 0), 0.95)         # per-tree annual mortality
}

# Drop-in for HiGy.R calc_mortality(): compute stand SDI from the plot summary
#   (sdi = tph.plot * (qmd/25)^1.6), then either
#   (a) uniform stand rate:
#       surv = koa.SURV.calibrated(sdi = sdi, planted = stand$planted)
#       dexpf = expf * (1 - surv) * mort.mult
#   (b) size-allocated (recommended for individual-tree realism):
#       m_stand = 1 - koa.SURV.calibrated(sdi = sdi, planted = stand$planted)
#       p_tree  = koa.SURV.allocate(dbh, expf, qmd, m_stand)
#       dexpf   = expf * p_tree * mort.mult
# If any calling code hard-codes SDImax = 500, change it to 1350 in the SAME
# commit, and make sure koa_ingrowth.R is on the same constant; the two files
# must not disagree about the density scale.
# Refine the constants as Kahikinui, KMR, and Kualoa remeasurements accrue.
