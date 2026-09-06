# ==============================================================================
# koa_prediction_functions.R
#
# Finalized prediction functions for the Acacia koa individual-tree
# growth and survival model system.
# Weiskittel, Sprecher, Gottesman & Rice (2025)
#
# All functions use population-average fixed-effects parameters derived from
# the best-fitting model variants selected during model development:
#
#   HT   -- Chapman-Richards with BYI and relative diameter (rDBH)
#   HCB  -- Logistic with slenderness, stand density, and BYI
#   dDBH -- Log-linear with BAL, CR, BAPH, Planted, BYI  (BYI variant)
#   dHT  -- Log-linear with BAL, CR, BAPH, Planted, BYI  (BYI variant)
#   SURV -- Complementary log-log GLM with HT, CR, rHT, BYI
#   BAL  -- Logistic function of relative diameter (rDBH = DBH/QMD)
#
# Units: DBH in cm, HT/HCB in m, BAPH/BAL in m2 ha-1,
#        TPH in trees ha-1, BYI in Mg ha-1, rain in mm, temp in deg C.
#        YIP (years in period) is used only in the multi-year wrapper functions.
#
# Smearing bias correction is applied to both increment functions:
#   CF_dDBH = 1.026  (exp(0.5 * s2) where s2 = WLS residual variance)
#   CF_dHT  = 1.030
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. STATIC HEIGHT MODEL
#    Form: HT = (a0 + a1*BYI/100) * (1 - exp(-b*DBH))^c
#               * exp(g1*log(BAPH+1) + g2*rDBH)
#    Data: AK_HT.csv, n = 10,060; R2 = 0.787; RMSE = 2.47 m
#    Note: rDBH = DBH/QMD; the negative g1 reflects stand-level density
#          suppression rather than individual suppression (population-average
#          formulation).
# ------------------------------------------------------------------------------

koa.HT <- function(DBH, BAPH, QMD, BYI = 264) {
  # DBH  : subject tree DBH (cm)
  # BAPH : stand basal area per hectare (m2 ha-1)
  # QMD  : stand quadratic mean diameter (cm)
  # BYI  : Biomass Yield Index (Mg ha-1); default = 264 (medium site)

  a0 <- 19.832
  a1 <-  0.106
  b  <-  0.044
  c  <-  0.863
  g1 <- -0.198
  g2 <-  0.479

  rDBH <- DBH / pmax(QMD, 0.1)
  HT   <- (a0 + a1 * BYI / 100) *
           (1 - exp(-b * DBH))^c *
           exp(g1 * log(BAPH + 1) + g2 * rDBH)
  return(pmax(HT, 1.37))   # minimum = breast height
}


# ------------------------------------------------------------------------------
# 2. HEIGHT TO CROWN BASE MODEL
#    Form: HCB = HT / (1 + exp(-eta))
#    eta  = b0 + b1*sqrt(HT/100) + b2*log(HT/DBH) + b3*sqrt(BAL*BAPH+1)
#               + b4*log(BAPH+1) + b5*log(BYI)
#    Data: AK_HCB.csv, n = 360; R2 = 0.763 (based on fitted vs observed CR)
#    Note: slenderness encoded as log(HT/DBH); BYI effect (b5 < 0) reflects
#          deeper crowns on more productive sites.
# ------------------------------------------------------------------------------

koa.HCB <- function(DBH, HT, BAPH, BAL, BYI = 264) {
  # DBH  : subject tree DBH (cm)
  # HT   : total height (m)
  # BAPH : stand basal area (m2 ha-1)
  # BAL  : basal area of trees larger than subject tree (m2 ha-1)
  # BYI  : Biomass Yield Index (Mg ha-1)

  b0 <-  1.162
  b1 <- -0.520
  b2 <- -0.257
  b3 <- -0.0147
  b4 <- -0.119
  b5 <- -0.221

  slender <- log(pmax(HT / pmax(DBH, 0.1), 0.5))
  eta     <- b0 + b1 * sqrt(pmax(HT / 100, 0)) + b2 * slender +
             b3 * sqrt(pmax(BAL * BAPH, 0) + 1) +
             b4 * log(BAPH + 1) + b5 * log(BYI)
  HCB     <- HT / (1 + exp(-eta))
  return(pmin(pmax(HCB, 0), 0.95 * HT))
}


# ------------------------------------------------------------------------------
# 3. BAL ALLOCATION MODEL
#    Form: BAL_frac(rDBH) = 1 / (1 + exp(-1.842 + 3.956 * rDBH))
#    where rDBH = DBH/QMD; BAL = BAPH * BAL_frac
#    Data: AK_TREE_incr.csv, n = 13,492 tree-period records
#    Interpretation:
#      rDBH = 1.0 (average tree) -> BAL_frac = 0.108
#      rDBH = 1.5 (dominant)     -> BAL_frac = 0.016
# ------------------------------------------------------------------------------

koa.BAL.fraction <- function(rDBH) {
  # rDBH : relative diameter = DBH / QMD
  return(1.0 / (1.0 + exp(-1.842 + 3.956 * rDBH)))
}

koa.BAL <- function(DBH, QMD, BAPH) {
  # Convenience wrapper: returns BAL (m2 ha-1) for a single tree.
  rDBH <- DBH / pmax(QMD, 0.1)
  return(BAPH * koa.BAL.fraction(rDBH))
}


# ------------------------------------------------------------------------------
# 4. ANNUAL DIAMETER INCREMENT
#    Form: dDBH_ann = exp(lp) * CF_dDBH
#    lp   = b0 + b1*log(DBH+1) + b2*DBH + b3*BAL^2/log(DBH+5)
#               + b4*log(BAL+1) + b5*log(CR) + b6*sqrt(BAPH*DBH)
#               + b7*Planted*DBH + b8*log(BYI)
#    Data: dDBH.csv, n = 6,209; WLS weights = 1/sqrt(YIP)
#    Performance: RMSE = 1.424 cm yr-1, R2 = 0.298
#    Note: Multi-year prediction -- use koa.dDBH.period() wrapper below.
# ------------------------------------------------------------------------------

CF_dDBH <- 1.026   # Duan smearing correction factor

koa.dDBH.annual <- function(DBH, BAL, CR, BAPH, Planted = 0, BYI = 264) {
  # DBH     : start-of-year DBH (cm)
  # BAL     : basal area of larger trees at start of year (m2 ha-1)
  # CR      : crown ratio (live crown length / total height), 0-1
  # BAPH    : stand basal area (m2 ha-1)
  # Planted : 1 if planted stand, 0 if natural
  # BYI     : Biomass Yield Index (Mg ha-1)

  b0 <- -2.4704737
  b1 <-  0.2072221
  b2 <- -0.0159616
  b3 <- -0.0016893
  b4 <- -0.2972574
  b5 <- -0.4470330
  b6 <- -0.0158403
  b7 <-  0.0188938
  b8 <-  0.4530166

  lp <- b0 + b1 * log(DBH + 1) + b2 * DBH +
        b3 * BAL^2 / log(DBH + 5) + b4 * log(BAL + 1) +
        b5 * log(pmax(CR, 0.01)) + b6 * sqrt(BAPH * DBH) +
        b7 * Planted * DBH + b8 * log(BYI)
  return(pmin(exp(lp) * CF_dDBH, 6))   # upper clip at 6 cm yr-1
}


koa.dDBH.period <- function(DBH.0, BAL.0, BAL.1, CR.0, CR.1,
                             BAPH.0, BAPH.1, Planted = 0, BYI = 264, YIP) {
  # Multi-year diameter growth by annual stepping with linear interpolation
  # of competition and crown variables across the period.
  # Returns predicted total diameter growth (cm) over YIP years.
  #
  # DBH.0, BAL.0, BAL.1, CR.0, CR.1, BAPH.0, BAPH.1 :
  #   Start- and end-of-period values (vectors OK).
  # YIP : years in period (vector OK, may vary by observation).

  n    <- length(DBH.0)
  dDBH <- numeric(n)

  for (i in seq_len(n)) {
    d      <- DBH.0[i]
    bal.c  <- BAL.0[i];  cr.c  <- CR.0[i];  bapa.c <- BAPH.0[i]
    bal.gr <- (BAL.1[i]  - BAL.0[i])  / YIP[i]
    cr.gr  <- (CR.1[i]   - CR.0[i])   / YIP[i]
    bapa.gr<- (BAPH.1[i] - BAPH.0[i]) / YIP[i]
    for (t in seq_len(YIP[i])) {
      gr    <- koa.dDBH.annual(d, bal.c, cr.c, bapa.c, Planted[i], BYI[i])
      d     <- d + gr
      bal.c <- bal.c  + bal.gr
      cr.c  <- cr.c   + cr.gr
      bapa.c<- bapa.c + bapa.gr
    }
    dDBH[i] <- d - DBH.0[i]
  }
  return(dDBH)
}


# ------------------------------------------------------------------------------
# 5. ANNUAL HEIGHT INCREMENT
#    Form: dHT_ann = exp(lp) * CF_dHT
#    lp   = b0 + b1*log(HT+1) + b2*HT + b3*BAL^2/log(HT+5)
#               + b4*log(BAL+1) + b5*log(CR) + b6*sqrt(BAPH*HT)
#               + b7*sqrt(Planted*HT) + b8*log(BYI)
#    Data: dHT.csv, n = 5,012; WLS weights = 1/sqrt(YIP)
#    Performance: RMSE = 1.082 m yr-1, R2 = 0.218
#    Note: Multi-year prediction -- use koa.dHT.period() wrapper below.
# ------------------------------------------------------------------------------

CF_dHT <- 1.030   # Duan smearing correction factor

koa.dHT.annual <- function(HT, BAL, CR, BAPH, Planted = 0, BYI = 264) {
  # HT      : start-of-year total height (m)
  # BAL     : basal area of larger trees (m2 ha-1)
  # CR      : crown ratio
  # BAPH    : stand basal area (m2 ha-1)
  # Planted : 1 if planted stand, 0 if natural
  # BYI     : Biomass Yield Index (Mg ha-1)

  b0 <- -3.382162
  b1 <-  0.272454
  b2 <- -0.105319
  b3 <- -0.000829
  b4 <- -0.071718
  b5 <- -1.483889
  b6 <-  0.033035
  b7 <-  0.017887
  b8 <-  0.433224

  lp <- b0 + b1 * log(HT + 1) + b2 * HT +
        b3 * BAL^2 / log(HT + 5) + b4 * log(BAL + 1) +
        b5 * log(pmax(CR, 0.01)) + b6 * sqrt(BAPH * HT) +
        b7 * sqrt(Planted * HT) + b8 * log(BYI)
  return(pmin(exp(lp) * CF_dHT, 4))   # upper clip at 4 m yr-1
}


koa.dHT.period <- function(HT.0, BAL.0, BAL.1, CR.0, CR.1,
                            BAPH.0, BAPH.1, Planted = 0, BYI = 264, YIP) {
  # Multi-year height growth by annual stepping with linear interpolation.
  # Returns predicted total height growth (m) over YIP years.

  n   <- length(HT.0)
  dHT <- numeric(n)

  for (i in seq_len(n)) {
    h      <- HT.0[i]
    bal.c  <- BAL.0[i];  cr.c  <- CR.0[i];  bapa.c <- BAPH.0[i]
    bal.gr <- (BAL.1[i]  - BAL.0[i])  / YIP[i]
    cr.gr  <- (CR.1[i]   - CR.0[i])   / YIP[i]
    bapa.gr<- (BAPH.1[i] - BAPH.0[i]) / YIP[i]
    for (t in seq_len(YIP[i])) {
      gr    <- koa.dHT.annual(h, bal.c, cr.c, bapa.c, Planted[i], BYI[i])
      h     <- h + gr
      bal.c <- bal.c  + bal.gr
      cr.c  <- cr.c   + cr.gr
      bapa.c<- bapa.c + bapa.gr
    }
    dHT[i] <- h - HT.0[i]
  }
  return(dHT)
}


# ------------------------------------------------------------------------------
# 6. TREE SURVIVAL MODEL
#    Form: P(alive over YIP years) = exp(-exp(lp) * YIP)
#          where lp = eta below (complementary log-log, population-average GLM)
#    eta  = b0 + b1*HT + b2*log(HT) + b3*rHT + b4*log(CR)
#               + b5*log(HT/DBH) + b6*log(BYI/100) + b7*(BYI/100)/1000
#    Data: AK_SURV.csv, n = 6,489, Deaths = 162 (2.5%)
#    Performance: AUC = 0.938; CV AUC = 0.931 (SD = 0.021)
#    Reference: Weiskittel et al. (2025), fitted as GLM with cloglog link
#               and log(YIP) offset. BYI variant (model S-B0, 8 parameters).
# ------------------------------------------------------------------------------

koa.surv <- function(DBH, HT, CR, rHT, BYI = 264, YIP = 1) {
  # DBH : start-of-period DBH (cm); used for height:diameter ratio only
  # HT  : total height (m)
  # CR  : crown ratio
  # rHT : relative height = HT / maximum stand height (set to 0.50 for
  #       average-tree stand-level simulation)
  # BYI : Biomass Yield Index (Mg ha-1)
  # YIP : years in period; returns period-level survival probability

  b0 <- 18.133
  b1 <-  0.199
  b2 <- -5.718
  b3 <-  7.640
  b4 <- 15.678
  b5 <- -3.396
  b6 <-  3.039
  b7 <- -25.102

  BYI.s <- BYI / 100
  HD    <- HT / pmax(DBH / 100, 0.01)   # height:diameter ratio (m m-1)

  eta   <- b0 + b1 * HT + b2 * log(pmax(HT, 0.5)) +
           b3 * rHT + b4 * log(pmax(CR, 0.01)) +
           b5 * log(pmax(HD, 1)) + b6 * log(BYI.s) + b7 * (BYI.s / 1000)

  # Period-level survival: P = exp(-exp(eta) * YIP)
  p.surv <- exp(-exp(eta) * YIP)
  return(pmin(pmax(p.surv, 0), 1))
}


# ==============================================================================
# CONVENIENCE WRAPPER: stand-table projection for a single cohort
# ==============================================================================

koa.project <- function(BYI = 264, Planted = 0,
                         init.DBH = NULL, init.TPH = NULL,
                         init.age = 5,   max.age  = 110,
                         SDI.max  = 500) {
  # Single-cohort stand projection using all sub-models above.
  # Returns a data frame with one row per year.
  #
  # BYI      : Biomass Yield Index (Mg ha-1)
  # Planted  : 1 = planted, 0 = natural
  # init.DBH : starting QMD (cm); defaults depend on origin
  # init.TPH : starting density (trees ha-1)
  # SDI.max  : maximum SDI for self-thinning threshold

  DBH.MAX  <- if (Planted) 60 else 90
  if (is.null(init.DBH)) init.DBH <- if (Planted) 3.5 else 2.5
  if (is.null(init.TPH)) init.TPH <- if (Planted) 1600 else 1000

  ages <- seq(init.age, max.age)
  n    <- length(ages)
  res  <- data.frame(
    Age  = ages,
    QMD  = NA_real_, HT  = NA_real_, HTd = NA_real_,
    BAPH = NA_real_, TPH = NA_real_, CR  = NA_real_,
    SDI  = NA_real_, HCB = NA_real_, VOL = NA_real_
  )

  DBH  <- init.DBH
  TPH  <- init.TPH
  BAPH <- TPH * pi / 4 * (DBH / 100)^2
  QMD  <- DBH
  HT   <- koa.HT(DBH, pmax(BAPH, 0.1), QMD, BYI)
  CR   <- 0.65

  for (s in seq_len(n)) {
    SDI     <- TPH * (QMD / 25.4)^1.605
    BAL.avg <- BAPH * koa.BAL.fraction(1.0)
    BAL.dom <- BAPH * koa.BAL.fraction(1.5)
    HCB     <- koa.HCB(DBH, HT, BAPH, BAL.avg, BYI)
    CR      <- if (HT > 0.1) pmax(0.20, (HT - HCB) / HT) else 0.65
    VOL     <- BAPH * HT * 0.40   # form factor = 0.40

    res[s, "QMD"]  <- QMD
    res[s, "HT"]   <- HT
    res[s, "HTd"]  <- HT * 1.17   # approximate H40
    res[s, "BAPH"] <- BAPH
    res[s, "TPH"]  <- TPH
    res[s, "CR"]   <- CR
    res[s, "SDI"]  <- SDI
    res[s, "HCB"]  <- HCB
    res[s, "VOL"]  <- VOL

    if (s == n) break

    # Self-thinning
    SDI.pct <- SDI / SDI.max
    if (SDI.pct > 0.60) {
      target  <- 0.55 / SDI.pct
      TPH.new <- TPH * target
      size.adj<- (TPH / pmax(TPH.new, 1))^0.08
      DBH  <- DBH * size.adj
      QMD  <- DBH
      TPH  <- TPH.new
      BAPH <- TPH * pi / 4 * (QMD / 100)^2
    }

    # Survival (annual)
    p.ann <- koa.surv(DBH, HT, CR, rHT = 0.50, BYI = BYI, YIP = 1)
    TPH   <- TPH * p.ann
    if (TPH < 5) {
      res[(s + 1):n, ] <- res[s, ]
      res[(s + 1):n, "Age"] <- ages[(s + 1):n]
      break
    }
    BAPH <- TPH * pi / 4 * (QMD / 100)^2

    # Growth (annual)
    dDBH <- koa.dDBH.annual(DBH, BAL.avg, CR, BAPH, Planted, BYI)
    dHT  <- koa.dHT.annual( HT,  BAL.avg, CR, BAPH, Planted, BYI)
    DBH  <- DBH + dDBH
    QMD  <- DBH
    HT   <- HT  + dHT
    # Blend dynamic HT with static H-D expectation
    HT.static <- koa.HT(DBH, BAPH, QMD, BYI)
    HT   <- 0.65 * HT + 0.35 * HT.static
    HT   <- min(HT, 25 + BYI / 55)   # biological ceiling
    DBH  <- min(DBH, DBH.MAX)
    QMD  <- DBH
    BAPH <- TPH * pi / 4 * (QMD / 100)^2
  }

  return(res)
}


# ==============================================================================
# EXAMPLE USAGE
# ==============================================================================

if (FALSE) {

  # Project a medium-quality natural stand
  nat.med <- koa.project(BYI = 264, Planted = 0)
  head(nat.med)

  # Project three site classes and plot BAPH trajectories
  byi.vals <- c(100, 264, 450)
  cols      <- c("#2166ac", "#1a9850", "#d73027")
  plot(NA, xlim = c(5, 110), ylim = c(0, 30),
       xlab = "Stand age (yr)", ylab = "Basal area (m2 ha-1)")
  for (i in seq_along(byi.vals)) {
    df <- koa.project(BYI = byi.vals[i], Planted = 0)
    lines(df$Age, df$BAPH, col = cols[i], lwd = 2)
  }
  legend("topleft", legend = paste("BYI =", byi.vals),
         col = cols, lwd = 2, bty = "n")

  # Predict height and HCB for individual trees
  koa.HT(DBH = 30, BAPH = 20, QMD = 25, BYI = 264)
  koa.HCB(DBH = 30, HT = 18, BAPH = 20, BAL = 5, BYI = 264)

  # Predict annual growth for an individual tree
  koa.dDBH.annual(DBH = 20, BAL = 5, CR = 0.55, BAPH = 15, Planted = 0, BYI = 264)
  koa.dHT.annual( HT  = 14, BAL = 5, CR = 0.55, BAPH = 15, Planted = 0, BYI = 264)

  # Predict 5-year survival probability
  koa.surv(DBH = 20, HT = 14, CR = 0.55, rHT = 0.50, BYI = 264, YIP = 5)

}
