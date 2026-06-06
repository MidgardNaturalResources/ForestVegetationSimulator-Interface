# tune_survival.R  (base R)
# Data-driven calibration of koa annual mortality. Uses an exposure-based
# (actuarial) rate: annual mortality = deaths / tree-years (sum YIP), and a
# Poisson hazard-rate GLM with log(YIP) offset (stable for rare events, gives
# proper annual rates and interpretable rate ratios). Focus: origin, BYI, density.
d <- read.csv("/sessions/elegant-lucid-cannon/mnt/Koa/KoaDatasets_12152025/AK.SURV.csv")
d$died <- as.integer(d$Status.1 != "live")
keep <- with(d, DBH.0>0 & Status.0=="live" & YIP>0 & is.finite(DBH.0) & is.finite(BYI) &
               is.finite(BAPH.0) & is.finite(SDI.0))
d <- d[keep,]
d$exposure <- pmax(d$YIP, 0.1)
d$RD <- d$SDI.0 / 500          # relative density (SDImax ~ 500)
rate <- function(sub) c(n=nrow(sub), deaths=sum(sub$died),
                        tree_yrs=round(sum(sub$exposure),0),
                        ann_mort=round(sum(sub$died)/sum(sub$exposure),4))
cat("OVERALL annual mortality (deaths/tree-years):\n"); print(rate(d))
cat("\nBY ORIGIN (0=nat,1=plt):\n")
for(o in c(0,1)){cat(" origin",o,": "); print(rate(d[d$Planted==o,]))}

cat("\nBY DBH CLASS:\n")
d$dclass <- cut(d$DBH.0, c(0,5,10,20,40,200), labels=c("0-5","5-10","10-20","20-40","40+"))
for(l in levels(d$dclass)){s<-d[d$dclass==l,]; if(nrow(s)>20){cat(sprintf(" DBH %-6s",l)); print(rate(s))}}

cat("\nBY BYI TERCILE:\n")
d$byiclass <- cut(d$BYI, quantile(d$BYI,c(0,.33,.67,1),na.rm=T), include.lowest=T)
for(l in levels(d$byiclass)){s<-d[d$byiclass==l,]; cat(sprintf(" BYI %-14s",l)); print(rate(s))}

cat("\nBY RELATIVE DENSITY (SDI/500):\n")
d$rdclass <- cut(d$RD, c(0,.3,.55,.8,3), labels=c("<.3",".3-.55",".55-.8",">.8"))
for(l in levels(d$rdclass)){s<-d[!is.na(d$rdclass)&d$rdclass==l,]; if(nrow(s)>20){cat(sprintf(" RD %-7s",l)); print(rate(s))}}

cat("\nBY ORIGIN x BYI tercile (annual mort):\n")
print(round(tapply(d$died,list(origin=d$Planted,byi=d$byiclass),sum)/
            tapply(d$exposure,list(origin=d$Planted,byi=d$byiclass),sum),4))

# ---- Poisson hazard-rate model (parsimonious, stable) ----
cat("\n=== Poisson rate model: died ~ log(DBH) + Planted + RD + log(BYI), offset log(YIP) ===\n")
m <- glm(died ~ I(log(DBH.0)) + Planted + RD + I(log(BYI)),
         family=poisson(link="log"), offset=log(exposure), data=d)
print(round(summary(m)$coefficients,4))
cat("\nRate ratios (exp coef):\n"); print(round(exp(coef(m)),4))
# predicted annual rate range
d$pred_ann <- exp(predict(m, newdata=transform(d, exposure=1)) )  # offset log(1)=0
cat("\npredicted annual mortality: median",round(median(d$pred_ann),4),
    " mean",round(mean(d$pred_ann),4)," 95th",round(quantile(d$pred_ann,.95),4),"\n")
cat("predicted by origin: nat",round(mean(d$pred_ann[d$Planted==0]),4),
    " plt",round(mean(d$pred_ann[d$Planted==1]),4),"\n")
# also without BYI (robustness given limited data)
m2 <- glm(died ~ I(log(DBH.0)) + Planted + RD, family=poisson, offset=log(exposure), data=d)
cat("\n=== Poisson WITHOUT BYI ===\n"); print(round(exp(coef(m2)),4))
cat("AIC with BYI:",round(AIC(m),1)," without BYI:",round(AIC(m2),1),"\n")

# write coefficients
co <- coef(m); con<-file("/sessions/elegant-lucid-cannon/mnt/outputs/koa_harness/results/survival_poisson_coef.csv","w")
cat("term,estimate\n",file=con); for(t in names(co)) cat(sprintf("%s,%.8f\n",t,co[t]),file=con); close(con)
co2 <- coef(m2); con<-file("/sessions/elegant-lucid-cannon/mnt/outputs/koa_harness/results/survival_poisson_nobyi_coef.csv","w")
cat("term,estimate\n",file=con); for(t in names(co2)) cat(sprintf("%s,%.8f\n",t,co2[t]),file=con); close(con)
cat("\nWrote survival_poisson_coef.csv (+ no-BYI)\n")
