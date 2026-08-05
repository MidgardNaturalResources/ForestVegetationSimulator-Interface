# check_origin_byi.R -- does survival change by origin and/or BYI once density
# (self-thinning) is accounted for? Exposure-based annual mortality rates.
d <- read.csv("/sessions/elegant-lucid-cannon/mnt/Koa/KoaDatasets_12152025/AK.SURV.csv")
d$died <- as.integer(d$Status.1 != "live")
d <- d[with(d, DBH.0>0 & Status.0=="live" & YIP>0 & is.finite(BYI) & is.finite(SDI.0)),]
d$exp <- pmax(d$YIP,0.1); d$RD <- d$SDI.0/500
rate <- function(s) if(nrow(s)<15) c(n=nrow(s),deaths=sum(s$died),ann=NA) else
  c(n=nrow(s),deaths=sum(s$died),tree_yrs=round(sum(s$exp)),ann=round(sum(s$died)/sum(s$exp),4))

cat("=== BELOW self-thinning onset (RD < 0.65): background mortality ===\n")
lo <- d[d$RD<0.65,]
cat("overall: "); print(rate(lo))
cat("by origin (0=nat,1=plt):\n"); for(o in c(0,1)){cat(" o",o,": ");print(rate(lo[lo$Planted==o,]))}
cat("by BYI tercile (within low density):\n")
lo$bt <- cut(lo$BYI, quantile(lo$BYI,c(0,.33,.67,1),na.rm=T), include.lowest=T)
for(l in levels(lo$bt)){cat(sprintf(" BYI %-14s",l)); print(rate(lo[lo$bt==l,]))}
cat("natural-only, by BYI tercile (low density):\n")
ln<-lo[lo$Planted==0,]; for(l in levels(lo$bt)){cat(sprintf(" nat BYI %-14s",l)); print(rate(ln[ln$bt==l,]))}

cat("\n=== ABOVE onset (RD >= 0.65): self-thinning zone, by origin ===\n")
hi <- d[d$RD>=0.65,]
cat("overall: "); print(rate(hi))
for(o in c(0,1)){cat(" o",o,": ");print(rate(hi[hi$Planted==o,]))}
cat("by RD band:\n")
hi$rb<-cut(hi$RD,c(.65,.8,1,3),labels=c(".65-.8",".8-1",">1"))
for(l in levels(hi$rb)){cat(sprintf(" RD %-7s",l));print(rate(hi[!is.na(hi$rb)&hi$rb==l,]))}

cat("\n=== Formal tests (Poisson rate, offset log exposure) ===\n")
cat("Background subset (RD<0.65): add origin? add BYI?\n")
m0<-glm(died~1,poisson,offset=log(exp),data=lo)
mO<-glm(died~Planted,poisson,offset=log(exp),data=lo)
mOB<-glm(died~Planted+I(log(BYI)),poisson,offset=log(exp),data=lo)
cat(sprintf(" null AIC %.1f | +origin AIC %.1f (dAIC %.1f) | +origin+BYI AIC %.1f (dAIC %.1f)\n",
  AIC(m0),AIC(mO),AIC(m0)-AIC(mO),AIC(mOB),AIC(mO)-AIC(mOB)))
cat(" origin coef (rate ratio): ",round(exp(coef(mO)["Planted"]),3),
    " | BYI logRR (given origin): ",round(coef(mOB)["I(log(BYI))"],3),
    " p=",signif(summary(mOB)$coef["I(log(BYI))","Pr(>|z|)"],3),"\n")
