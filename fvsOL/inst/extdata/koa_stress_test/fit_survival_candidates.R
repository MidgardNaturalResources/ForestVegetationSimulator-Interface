# fit_survival_candidates.R  (base R only)
# Fit and diagnose survival candidates on the current AK.SURV.csv, focused on
# operational stability (esp. plantations), not just AUC. Outputs coefficients
# to JSON-like CSV + diagnostics.
d <- read.csv("/sessions/elegant-lucid-cannon/mnt/Koa/KoaDatasets_12152025/AK.SURV.csv")
d$alive  <- as.integer(d$Status.1 == "live")
d$logYIP <- log(pmax(d$YIP, 0.1))
keep <- with(d, DBH.0>0 & Status.0=="live" & YIP>0 & is.finite(HT.0) & is.finite(CR.0) &
               is.finite(rHT.0) & is.finite(BYI) & CR.0>0 & HT.0>0)
d <- d[keep,]
cat("n =", nrow(d), " deaths =", sum(d$alive==0), " mort =", round(mean(1-d$alive),4), "\n")
cat("mortality by origin (0=nat,1=plt): "); print(round(tapply(1-d$alive, d$Planted, mean),4))

auc <- function(p,y){ r<-rank(p); (sum(r[y==1])-sum(y==1)*(sum(y==1)+1)/2)/(sum(y==1)*sum(y==0)) }

specs <- list(
  S1_published = list(link="cloglog", f=alive ~ HT.0 + I(log(HT.0)) + rHT.0 + I(log(CR.0)) +
                        I(log(HT.0/(DBH.0/100))) + I(log(BYI/100)) + I(BYI/1000)),
  S4_dbh_planted = list(link="cloglog", f=alive ~ DBH.0 + I(log(DBH.0)) +
                        I((BAL.0+1)/log(DBH.0+1)) + I(log(CR.0)) + I(sqrt(BAPH.0*DBH.0)) +
                        I(log(HT.0/(DBH.0/100))) + I(Planted*DBH.0)),
  S3_dbh_planted_byi = list(link="cloglog", f=alive ~ DBH.0 + I(log(DBH.0)) +
                        I((BAL.0+1)/log(DBH.0+1)) + I(log(CR.0)) + I(sqrt(BAPH.0*DBH.0)) +
                        I(log(HT.0/(DBH.0/100))) + I(Planted*DBH.0) + log(BYI) + I(BYI/1000)),
  Snew_logit = list(link="logit", f=alive ~ I(log(DBH.0)) + DBH.0 +
                        I((BAL.0+1)/log(DBH.0+1)) + I(log(CR.0)) + Planted + I(BYI/1000)),
  Snew2_logit = list(link="logit", f=alive ~ I(log(DBH.0)) +
                        I(log(pmax(BAL.0,0)+1)) + I(log(CR.0)) + Planted + I(log(BYI/100))),
  # cloglog (proper annualization) with DBH main effect + Planted + tempered terms
  Snew_cloglog = list(link="cloglog", f=alive ~ I(log(DBH.0)) + DBH.0 +
                        I(log(pmax(BAL.0,0)+1)) + I(log(CR.0)) + Planted + I(BYI/1000)),
  Snew2_cloglog = list(link="cloglog", f=alive ~ I(log(DBH.0)) +
                        I(log(pmax(BAL.0,0)+1)) + I(log(CR.0)) + Planted + I(log(BYI/100))),
  # MINIMAL, monotonic, stable: size + competition + origin (proper cloglog annualization)
  Smin_cloglog  = list(link="cloglog", f=alive ~ I(log(DBH.0)) +
                        I((BAL.0+1)/log(DBH.0+1)) + Planted),
  Smin2_cloglog = list(link="cloglog", f=alive ~ I(log(DBH.0)) + Planted),
  Smin3_cloglog = list(link="cloglog", f=alive ~ I(log(DBH.0)) +
                        I(log(pmax(BAL.0,0)+1)) + Planted + I(log(CR.0)))
)

res <- list()
for (nm in names(specs)) {
  sp <- specs[[nm]]
  fit <- suppressWarnings(glm(sp$f, data=d, family=binomial(link=sp$link), offset=d$logYIP))
  p <- suppressWarnings(predict(fit, type="response"))   # P(alive) at obs YIP
  # annual prob (YIP=1): refit linpred minus offset
  eta <- predict(fit, type="link") - d$logYIP
  pann <- if (sp$link=="cloglog") 1-exp(-exp(eta)) else plogis(eta)
  mort_ann <- 1 - pann
  res[[nm]] <- list(aic=AIC(fit), auc=auc(p,d$alive),
                    coef=coef(fit),
                    mort_ann_nat=mean(mort_ann[d$Planted==0]),
                    mort_ann_plt=mean(mort_ann[d$Planted==1]),
                    pct_surv_lt_080=mean(pann<0.80),
                    pct_surv_lt_050=mean(pann<0.50))
  cat(sprintf("\n%-20s link=%-8s AIC=%.0f AUC=%.3f | annual mort mean nat=%.3f plt=%.3f | median nat=%.3f plt=%.3f | %%pann<0.5=%.1f\n",
      nm, sp$link, AIC(fit), auc(p,d$alive), res[[nm]]$mort_ann_nat, res[[nm]]$mort_ann_plt,
      median(mort_ann[d$Planted==0]), median(mort_ann[d$Planted==1]), 100*res[[nm]]$pct_surv_lt_050))
}
cat("\nObserved annual mortality (approx): nat ~",
    round(mean((1-d$alive)[d$Planted==0])/mean(d$YIP[d$Planted==0]),4),
    " plt ~", round(mean((1-d$alive)[d$Planted==1])/mean(d$YIP[d$Planted==1]),4), "\n")

# write coefficients
sink("/sessions/elegant-lucid-cannon/mnt/outputs/koa_harness/results/survival_candidates_coef.txt")
for (nm in names(res)) { cat("=== ",nm," (link=",specs[[nm]]$link,") AUC=",round(res[[nm]]$auc,3)," ===\n",sep="")
  print(round(res[[nm]]$coef,6)); cat("\n") }
sink()
# machine-readable
con <- file("/sessions/elegant-lucid-cannon/mnt/outputs/koa_harness/results/survival_candidates.csv","w")
cat("variant,link,term,estimate\n", file=con)
for (nm in names(res)) for (t in names(res[[nm]]$coef))
  cat(sprintf("%s,%s,%s,%.8f\n", nm, specs[[nm]]$link, t, res[[nm]]$coef[t]), file=con)
close(con)
cat("\nWrote survival_candidates_coef.txt and survival_candidates.csv\n")
