################################################################################
# 2026.02.23
#  Koa (Acacia koa) Individual Tree Growth Model - Prediction Functions
#  Weiskittel, Sprecher, Gottesman & Rice (2025)
#  For use in growth and yield projections across Hawaiian koa forests
#
#  Component equations:
#    1. Total height (HT) - Chapman-Richards, with and without BYI
#    2. Height to crown base (HCB) - logistic, bounded by HT
#    3. Diameter increment (dDBH) - exponential, with and without BYI
#    4. Height increment (dHT) - exponential, with and without BYI
#    5. Annual survival probability (PS) - cloglog, with and without BYI
#
#  Units:  DBH (cm), HT (m), HCB (m), BAPH (m2/ha), BAL (m2/ha),
#          BYI (Mg/ha), dDBH (cm/yr), dHT (m/yr)
#
#  Variables:
#    DBH      - diameter at breast height (cm)
#    HT       - total height (m)
#    HCB      - height to crown base (m)
#    CR       - crown ratio = (HT - HCB) / HT
#    BAPH     - stand basal area per ha (m2/ha)
#    BAL      - basal area in larger trees (m2/ha)
#    QMD      - quadratic mean diameter (cm)
#    rDBH     - relative diameter = DBH / QMD
#    BYI      - Biomass Yield Index (Mg/ha); range ~52 to 618
#    Planted  - origin indicator: 1 = planted, 0 = natural
#    YIP      - years in projection interval (for survival)
################################################################################


# ==============================================================================
#  1. STATIC TOTAL HEIGHT EQUATION
#     HT = (a0 + a1 * BYI/100) * [1 - exp(-b * DBH)]^c *
#          exp(g1 * ln(BAPH+1) + g2 * rDBH)
# ==============================================================================

#' Predict total height (m) from DBH and stand attributes
#'
#' @param DBH   Diameter at breast height (cm)
#' @param BAPH  Stand basal area per ha (m2/ha)
#' @param BAL   Basal area in larger trees (m2/ha)
#' @param QMD   Quadratic mean diameter (cm)
#' @param BYI   Biomass Yield Index (Mg/ha). If NULL, uses basic model.
#' @return Predicted total height (m)
predict_HT <- function(DBH, BAPH, BAL, QMD, BYI = NULL) {

  # Relative diameter
  rDBH <- DBH / QMD

  if (is.null(BYI)) {
    # ---- Basic model (no site quality) ----
    # Parameters from nlme fit; n = 10,060; R2 = 0.793; RMSE = 2.41 m
    a0 <- 22.074     # SE = 0.521  asymptote intercept
    b  <-  0.041     # SE = 0.002  growth rate
    c  <-  0.851     # SE = 0.018  shape
    g1 <- -0.214     # SE = 0.031  BAPH competition
    g2 <-  0.487     # SE = 0.028  relative diameter

    HT <- a0 * (1 - exp(-b * DBH))^c *
          exp(g1 * log(BAPH + 1) + g2 * rDBH)
  } else {
    # ---- BYI-enhanced model ----
    # n = 10,060; R2 = 0.803; RMSE = 2.35 m; Bias = 0.01 m (0.06%)
    a0 <- 19.832     # SE = 0.614  asymptote intercept
    a1 <-  0.106     # SE = 0.018  BYI effect on asymptote
    b  <-  0.044     # SE = 0.002  growth rate
    c  <-  0.863     # SE = 0.019  shape
    g1 <- -0.198     # SE = 0.030  BAPH competition
    g2 <-  0.479     # SE = 0.027  relative diameter

    HT <- (a0 + a1 * BYI / 100) * (1 - exp(-b * DBH))^c *
          exp(g1 * log(BAPH + 1) + g2 * rDBH)
  }

  return(pmax(HT, 1.37))  # enforce minimum above breast height
}


# ==============================================================================
#  2. STATIC HEIGHT TO CROWN BASE EQUATION
#     HCB = HT * {1 + exp[b0 + b1*sqrt(HT/100) + b2*ln(HT/DBH)
#                         + b3*sqrt(BAL*BAPH+1) + b4*ln(BAPH+1)]}^-1
# ==============================================================================

#' Predict height to live crown base (m)
#'
#' @param DBH   Diameter at breast height (cm)
#' @param HT    Total height (m)
#' @param BAPH  Stand basal area per ha (m2/ha)
#' @param BAL   Basal area in larger trees (m2/ha)
#' @return Predicted HCB (m); crown ratio can be computed as (HT - HCB) / HT
predict_HCB <- function(DBH, HT, BAPH, BAL) {

  # Parameters from nlme fit; n = 359; R2 = 0.408; RMSE = 2.11 m
  b0 <- -1.233     # SE = 0.241  intercept
  b1 <-  0.874     # SE = 0.193  sqrt(HT/100)
  b2 <-  0.249     # SE = 0.054  ln(HT/DBH) slenderness
  b3 <-  0.0041    # SE = 0.0012 sqrt(BAL*BAPH+1) competition interaction
  b4 <-  0.342     # SE = 0.071  ln(BAPH+1) stand density

  eta <- b0 +
         b1 * sqrt(HT / 100) +
         b2 * log(HT / DBH) +
         b3 * sqrt(BAL * BAPH + 1) +
         b4 * log(BAPH + 1)

  HCB <- HT / (1 + exp(-eta))

  # Crown base cannot exceed 95% of height or fall below ground
  HCB <- pmin(HCB, 0.95 * HT)
  HCB <- pmax(HCB, 0)

  return(HCB)
}


# ==============================================================================
#  3. DIAMETER INCREMENT EQUATION (annual, cm/yr)
#     dDBH = exp{b0 + b1*ln(DBH+1) + b2*DBH + b3*BAL^2/ln(DBH+5)
#               + b4*ln(BAL+1) + b5*sqrt(BAPH*DBH) + b6*ln(BAL*BAPH+1)
#               [+ b7*ln(BYI) + b8*BYI/1000]
#               [+ b9*Planted]}
# ==============================================================================

#' Predict annual diameter increment (cm/yr)
#'
#' @param DBH     Diameter at breast height (cm)
#' @param BAPH    Stand basal area per ha (m2/ha)
#' @param BAL     Basal area in larger trees (m2/ha)
#' @param BYI     Biomass Yield Index (Mg/ha). If NULL, uses basic model.
#' @param Planted Origin indicator (1 = planted, 0 = natural)
#' @return Predicted annual diameter increment (cm/yr)
predict_dDBH <- function(DBH, BAPH, BAL, BYI = NULL, Planted = 0) {

  if (is.null(BYI)) {
    # ---- Basic model (no site quality) ----
    # n = 6,209; R2 = 0.224; RMSE = 1.53 cm/yr; Bias = -0.09 cm/yr (-6.6%); CV = 0.197
    b0 <- -0.814     # SE = 0.114  intercept
    b1 <-  0.623     # SE = 0.042  ln(DBH+1) size (positive limb)
    b2 <- -0.011     # SE = 0.001  DBH (negative limb, creates hump)
    b3 <- -0.00019   # SE = 0.000023 BAL^2/ln(DBH+5) asymmetric competition
    b4 <- -0.419     # SE = 0.038  ln(BAL+1) competition
    b5 <- -0.018     # SE = 0.003  sqrt(BAPH*DBH) density-size interaction
    b6 <-  0.062     # SE = 0.012  ln(BAL*BAPH+1) density amelioration
    b9 <-  0.291     # SE = 0.031  Planted origin additive effect

    lp <- b0 +
          b1 * log(DBH + 1) +
          b2 * DBH +
          b3 * BAL^2 / log(DBH + 5) +
          b4 * log(BAL + 1) +
          b5 * sqrt(BAPH * DBH) +
          b6 * log(BAL * BAPH + 1) +
          b9 * Planted

  } else {
    # ---- BYI-enhanced model ----
    # n = 6,209; R2 = 0.268; RMSE = 1.47 cm/yr; Bias = -0.08 cm/yr (-5.9%); CV = 0.241
    b0 <- -4.217     # SE = 0.381  intercept
    b1 <-  0.611     # SE = 0.041  ln(DBH+1)
    b2 <- -0.010     # SE = 0.001  DBH
    b3 <- -0.00018   # SE = 0.000022 BAL^2/ln(DBH+5)
    b4 <- -0.401     # SE = 0.037  ln(BAL+1)
    b5 <- -0.017     # SE = 0.003  sqrt(BAPH*DBH)
    b6 <-  0.059     # SE = 0.011  ln(BAL*BAPH+1)
    b7 <-  0.952     # SE = 0.089  ln(BYI) - positive BYI effect
    b8 <- -0.00112   # SE = 0.00021 BYI/1000 - curvilinear (peak ~420 Mg/ha)
    b9 <-  0.287     # SE = 0.030  Planted

    lp <- b0 +
          b1 * log(DBH + 1) +
          b2 * DBH +
          b3 * BAL^2 / log(DBH + 5) +
          b4 * log(BAL + 1) +
          b5 * sqrt(BAPH * DBH) +
          b6 * log(BAL * BAPH + 1) +
          b7 * log(BYI) +
          b8 * BYI / 1000 +
          b9 * Planted
  }

  dDBH <- exp(lp)

  # Biological bounds: 0 to 8 cm/yr (values >5 cm/yr are very rare)
  dDBH <- pmin(pmax(dDBH, 0), 8)

  return(dDBH)
}


# ==============================================================================
#  4. HEIGHT INCREMENT EQUATION (annual, m/yr)
#     dHT = exp{b0 + b1*ln(HT+1) + b2*HT + b3*BAL^2/ln(DBH+5)
#              + b4*ln(BAL+1) + b5*sqrt(BAPH*DBH) + b6*ln(BAL*BAPH+1)
#              [+ b7*ln(BYI) + b8*BYI/1000]
#              [+ b9*Planted]}
# ==============================================================================

#' Predict annual height increment (m/yr)
#'
#' @param DBH     Diameter at breast height (cm)
#' @param HT      Total height (m)
#' @param BAPH    Stand basal area per ha (m2/ha)
#' @param BAL     Basal area in larger trees (m2/ha)
#' @param BYI     Biomass Yield Index (Mg/ha). If NULL, uses basic model.
#' @param Planted Origin indicator (1 = planted, 0 = natural)
#' @return Predicted annual height increment (m/yr)
predict_dHT <- function(DBH, HT, BAPH, BAL, BYI = NULL, Planted = 0) {

  if (is.null(BYI)) {
    # ---- Basic model (no site quality) ----
    # n = 5,012; R2 = 0.176; RMSE = 1.15 m/yr; Bias = -0.03 m/yr (-2.1%); CV = 0.163
    b0 <- -0.427     # SE = 0.097  intercept
    b1 <-  1.087     # SE = 0.071  ln(HT+1) size (positive limb)
    b2 <- -0.042     # SE = 0.004  HT (negative limb, creates hump)
    b3 <- -0.00011   # SE = 0.000018 BAL^2/ln(DBH+5)
    b4 <- -0.311     # SE = 0.034  ln(BAL+1)
    b5 <- -0.016     # SE = 0.003  sqrt(BAPH*DBH)
    b6 <-  0.047     # SE = 0.010  ln(BAL*BAPH+1)
    b9 <-  0.348     # SE = 0.035  Planted

    lp <- b0 +
          b1 * log(HT + 1) +
          b2 * HT +
          b3 * BAL^2 / log(DBH + 5) +
          b4 * log(BAL + 1) +
          b5 * sqrt(BAPH * DBH) +
          b6 * log(BAL * BAPH + 1) +
          b9 * Planted

  } else {
    # ---- BYI-enhanced model ----
    # n = 5,012; R2 = 0.213; RMSE = 1.09 m/yr; Bias = -0.02 m/yr (-1.4%); CV = 0.198
    b0 <- -3.891     # SE = 0.334  intercept
    b1 <-  1.074     # SE = 0.069  ln(HT+1)
    b2 <- -0.040     # SE = 0.004  HT
    b3 <- -0.00010   # SE = 0.000017 BAL^2/ln(DBH+5)
    b4 <- -0.298     # SE = 0.033  ln(BAL+1)
    b5 <- -0.015     # SE = 0.003  sqrt(BAPH*DBH)
    b6 <-  0.044     # SE = 0.010  ln(BAL*BAPH+1)
    b7 <-  0.869     # SE = 0.081  ln(BYI)
    b8 <- -0.00103   # SE = 0.00019 BYI/1000 (peak ~421 Mg/ha)
    b9 <-  0.343     # SE = 0.034  Planted

    lp <- b0 +
          b1 * log(HT + 1) +
          b2 * HT +
          b3 * BAL^2 / log(DBH + 5) +
          b4 * log(BAL + 1) +
          b5 * sqrt(BAPH * DBH) +
          b6 * log(BAL * BAPH + 1) +
          b7 * log(BYI) +
          b8 * BYI / 1000 +
          b9 * Planted
  }

  dHT <- exp(lp)

  # Biological bounds: 0 to 6 m/yr (values >4 m/yr are very rare)
  dHT <- pmin(pmax(dHT, 0), 6)

  return(dHT)
}


# ==============================================================================
#  5. ANNUAL SURVIVAL PROBABILITY (cloglog link)
#     PS = 1 - exp{-exp[b0 + b1*DBH + b2*ln(DBH^2) + b3*ln(HT/DBH+1)
#                       + b4*(BAL+1)/ln(DBH+1) + b5*ln(BAPH)
#                       [+ b6*BYI + b7*BYI^2]
#                       [+ b8*Planted]]}
#
#  For multi-year intervals: PS_interval = PS_annual ^ YIP
# ==============================================================================

#' Predict annual survival probability
#'
#' @param DBH     Diameter at breast height (cm)
#' @param HT      Total height (m)
#' @param BAPH    Stand basal area per ha (m2/ha)
#' @param BAL     Basal area in larger trees (m2/ha)
#' @param BYI     Biomass Yield Index (Mg/ha). If NULL, uses basic model.
#' @param Planted Origin indicator (1 = planted, 0 = natural)
#' @param YIP     Years in projection interval (default = 1 for annual)
#' @return Probability of survival over YIP years
predict_survival <- function(DBH, HT, BAPH, BAL, BYI = NULL,
                             Planted = 0, YIP = 1) {

  if (is.null(BYI)) {
    # ---- Basic cloglog model (no site quality) ----
    # n = 5,686; AUC = 0.876; Brier = 0.016; Youden's J = 0.829
    b0 <- -3.847     # SE = 0.418  intercept
    b1 <-  0.207     # SE = 0.031  DBH (size effect)
    b2 <- -0.686     # SE = 0.076  ln(DBH^2) (quadratic size)
    b3 <-  0.710     # SE = 0.099  ln(HT/DBH+1) slenderness
    b4 <- -0.024     # SE = 0.004  (BAL+1)/ln(DBH+1) asymmetric competition
    b5 <- -0.946     # SE = 0.127  ln(BAPH) stand density
    b8 <-  0.531     # SE = 0.089  Planted

    lp <- b0 +
          b1 * DBH +
          b2 * log(DBH^2) +
          b3 * log(HT / DBH + 1) +
          b4 * (BAL + 1) / log(DBH + 1) +
          b5 * log(BAPH) +
          b8 * Planted

  } else {
    # ---- BYI-enhanced cloglog model ----
    # n = 5,686; AUC = 0.893; Brier = 0.014
    b0 <- -2.914     # SE = 0.501  intercept
    b1 <-  0.198     # SE = 0.030  DBH
    b2 <- -0.661     # SE = 0.073  ln(DBH^2)
    b3 <-  0.683     # SE = 0.095  ln(HT/DBH+1)
    b4 <- -0.023     # SE = 0.004  (BAL+1)/ln(DBH+1)
    b5 <- -0.902     # SE = 0.122  ln(BAPH)
    b6 <-  0.00842   # SE = 0.00214 BYI (linear)
    b7 <- -0.0000185 # SE = 0.0000047 BYI^2 (quadratic; peak hazard ~227 Mg/ha)
    b8 <-  0.514     # SE = 0.086  Planted

    lp <- b0 +
          b1 * DBH +
          b2 * log(DBH^2) +
          b3 * log(HT / DBH + 1) +
          b4 * (BAL + 1) / log(DBH + 1) +
          b5 * log(BAPH) +
          b6 * BYI +
          b7 * BYI^2 +
          b8 * Planted
  }

  # cloglog: annual survival = exp(-exp(lp))
  PS_annual <- exp(-exp(lp))

  # Extend to multi-year interval
  PS_interval <- PS_annual^YIP

  PS_interval <- pmin(pmax(PS_interval, 0), 1)

  return(PS_interval)
}


# ==============================================================================
#  6. ITERATIVE ANNUAL PROJECTION FUNCTION
#     Projects a single tree forward one year, updating size and stand attributes
# ==============================================================================

#' Project a koa tree forward by one year
#'
#' @param tree  Named list or data frame row with: DBH, HT, BAPH, BAL, QMD,
#'              BYI (optional), Planted
#' @return Updated tree list with new DBH, HT, HCB, CR, and survival flag
project_tree_1yr <- function(tree) {

  BYI     <- tree[["BYI"]]      # NULL if not available
  Planted <- tree[["Planted"]]
  DBH     <- tree[["DBH"]]
  HT      <- tree[["HT"]]
  BAPH    <- tree[["BAPH"]]
  BAL     <- tree[["BAL"]]
  QMD     <- tree[["QMD"]]

  # Growth increments
  dDBH <- predict_dDBH(DBH, BAPH, BAL, BYI, Planted)
  dHT  <- predict_dHT(DBH, HT, BAPH, BAL, BYI, Planted)
  PS   <- predict_survival(DBH, HT, BAPH, BAL, BYI, Planted, YIP = 1)

  # Update size
  DBH_new <- DBH + dDBH
  HT_new  <- HT  + dHT
  HCB_new <- predict_HCB(DBH_new, HT_new, BAPH, BAL)
  CR_new  <- (HT_new - HCB_new) / HT_new

  # Stochastic survival (set alive = FALSE with probability 1 - PS)
  alive <- runif(1) < PS

  return(list(
    DBH     = DBH_new,
    HT      = HT_new,
    HCB     = HCB_new,
    CR      = CR_new,
    BAPH    = BAPH,    # stand-level update handled externally
    BAL     = BAL,
    QMD     = QMD,
    BYI     = BYI,
    Planted = Planted,
    alive   = alive,
    PS_annual = PS
  ))
}


# ==============================================================================
#  7. WORKED EXAMPLE
# ==============================================================================

cat("\n=== KOA GROWTH MODEL - WORKED PREDICTION EXAMPLES ===\n\n")

# --- Example tree attributes ---
ex_DBH     <- 25.0    # cm
ex_HT      <- 18.0    # m
ex_BAPH    <- 22.0    # m2/ha
ex_BAL     <- 12.0    # m2/ha
ex_QMD     <- 30.0    # cm
ex_BYI     <- 264.0   # Mg/ha (medium productivity)
ex_Planted <- 0       # natural origin

cat("--- Input tree attributes ---\n")
cat(sprintf("  DBH = %.1f cm,  HT = %.1f m\n", ex_DBH, ex_HT))
cat(sprintf("  BAPH = %.1f m2/ha,  BAL = %.1f m2/ha\n", ex_BAPH, ex_BAL))
cat(sprintf("  QMD = %.1f cm,  BYI = %.1f Mg/ha\n", ex_QMD, ex_BYI))
cat(sprintf("  Planted = %d (0 = natural)\n\n", ex_Planted))

# --- Static height predictions ---
HT_basic <- predict_HT(ex_DBH, ex_BAPH, ex_BAL, ex_QMD)
HT_byi   <- predict_HT(ex_DBH, ex_BAPH, ex_BAL, ex_QMD, BYI = ex_BYI)

cat("--- Total height predictions ---\n")
cat(sprintf("  Basic model:       %.2f m\n", HT_basic))
cat(sprintf("  BYI-enhanced:      %.2f m\n\n", HT_byi))

# --- Crown base ---
HCB_pred <- predict_HCB(ex_DBH, ex_HT, ex_BAPH, ex_BAL)
CR_pred  <- (ex_HT - HCB_pred) / ex_HT

cat("--- Crown base prediction ---\n")
cat(sprintf("  HCB = %.2f m,  Crown ratio = %.3f\n\n", HCB_pred, CR_pred))

# --- Diameter increment ---
dDBH_basic <- predict_dDBH(ex_DBH, ex_BAPH, ex_BAL)
dDBH_byi   <- predict_dDBH(ex_DBH, ex_BAPH, ex_BAL, BYI = ex_BYI)
dDBH_plt   <- predict_dDBH(ex_DBH, ex_BAPH, ex_BAL, BYI = ex_BYI, Planted = 1)

cat("--- Diameter increment predictions ---\n")
cat(sprintf("  Basic (natural):   %.3f cm/yr\n", dDBH_basic))
cat(sprintf("  BYI (natural):     %.3f cm/yr\n", dDBH_byi))
cat(sprintf("  BYI (planted):     %.3f cm/yr\n\n", dDBH_plt))

# --- Height increment ---
dHT_basic <- predict_dHT(ex_DBH, ex_HT, ex_BAPH, ex_BAL)
dHT_byi   <- predict_dHT(ex_DBH, ex_HT, ex_BAPH, ex_BAL, BYI = ex_BYI)

cat("--- Height increment predictions ---\n")
cat(sprintf("  Basic (natural):   %.3f m/yr\n", dHT_basic))
cat(sprintf("  BYI (natural):     %.3f m/yr\n\n", dHT_byi))

# --- Survival ---
PS_basic <- predict_survival(ex_DBH, ex_HT, ex_BAPH, ex_BAL)
PS_byi   <- predict_survival(ex_DBH, ex_HT, ex_BAPH, ex_BAL, BYI = ex_BYI)
PS_5yr   <- predict_survival(ex_DBH, ex_HT, ex_BAPH, ex_BAL, BYI = ex_BYI, YIP = 5)

cat("--- Survival probability predictions ---\n")
cat(sprintf("  Basic (1-yr):      %.4f\n", PS_basic))
cat(sprintf("  BYI (1-yr):        %.4f\n", PS_byi))
cat(sprintf("  BYI (5-yr):        %.4f\n\n", PS_5yr))

# --- Size-class summary table ---
cat("--- Predictions across DBH classes (BYI = 264, BAPH = 22, BAL = 8) ---\n")
cat(sprintf("%-8s %-8s %-8s %-8s %-8s %-8s\n",
            "DBH(cm)", "HT(m)", "HCB(m)", "CR", "dDBH", "PS(1yr)"))
cat(strrep("-", 56), "\n")

for (dbh in c(5, 10, 15, 20, 30, 40, 50, 60)) {
  ht   <- predict_HT(dbh, BAPH = 22, BAL = 8, QMD = 25, BYI = 264)
  hcb  <- predict_HCB(dbh, ht, BAPH = 22, BAL = 8)
  cr   <- (ht - hcb) / ht
  ddbh <- predict_dDBH(dbh, BAPH = 22, BAL = 8, BYI = 264)
  ps   <- predict_survival(dbh, ht, BAPH = 22, BAL = 8, BYI = 264)
  cat(sprintf("%-8.0f %-8.1f %-8.1f %-8.3f %-8.3f %-8.4f\n",
              dbh, ht, hcb, cr, ddbh, ps))
}


# ==============================================================================
#  8. SITE QUALITY (BYI) SENSITIVITY
# ==============================================================================

cat("\n--- Diameter increment across BYI classes (DBH = 25, BAL = 12, BAPH = 22) ---\n")
cat(sprintf("%-20s %-12s %-12s\n", "BYI class", "BYI (Mg/ha)", "dDBH (cm/yr)"))
cat(strrep("-", 48), "\n")

byi_classes <- list(
  c("Very Low",   42),
  c("Low",       120),
  c("Medium",    213),
  c("High",      347),
  c("Very High", 520)
)

for (cls in byi_classes) {
  byi_val <- as.numeric(cls[2])
  d <- predict_dDBH(25, BAPH = 22, BAL = 12, BYI = byi_val)
  cat(sprintf("%-20s %-12.0f %-12.3f\n", cls[1], byi_val, d))
}


# ==============================================================================
#  9. SIMPLE 10-YEAR PROJECTION EXAMPLE
# ==============================================================================

cat("\n--- 10-year deterministic projection (single tree) ---\n")
cat(sprintf("%-6s %-8s %-8s %-8s %-8s %-8s\n",
            "Year", "DBH(cm)", "HT(m)", "HCB(m)", "CR", "PS_ann"))
cat(strrep("-", 52), "\n")

# Initial conditions
cur_DBH  <- 15.0
cur_HT   <- 12.0
cur_BAPH <- 25.0
cur_BAL  <- 15.0
cur_QMD  <- 22.0
cur_BYI  <- 264.0
cur_Plt  <- 0

for (yr in 0:10) {
  cur_HCB <- predict_HCB(cur_DBH, cur_HT, cur_BAPH, cur_BAL)
  cur_CR  <- (cur_HT - cur_HCB) / cur_HT
  cur_PS  <- predict_survival(cur_DBH, cur_HT, cur_BAPH, cur_BAL, cur_BYI, cur_Plt)

  cat(sprintf("%-6d %-8.2f %-8.2f %-8.2f %-8.3f %-8.4f\n",
              yr, cur_DBH, cur_HT, cur_HCB, cur_CR, cur_PS))

  if (yr < 10) {
    # Update tree and stand (simplified: BAPH/BAL unchanged for illustration)
    cur_DBH <- cur_DBH + predict_dDBH(cur_DBH, cur_BAPH, cur_BAL, cur_BYI, cur_Plt)
    cur_HT  <- cur_HT  + predict_dHT(cur_DBH, cur_HT, cur_BAPH, cur_BAL, cur_BYI, cur_Plt)
  }
}

cat("\nNote: Stand-level variables (BAPH, BAL, QMD) held constant in this example.\n")
cat("In a full simulation, update BAPH, BAL, and QMD each year for all trees.\n")
cat("\nModel reference: Weiskittel, Sprecher, Gottesman & Rice (2025)\n")
cat("For questions contact: aaron.weiskittel@maine.edu\n")
