# check_interaction.R -- is there an estimable origin x BYI interaction?
d <- read.csv("/sessions/elegant-lucid-cannon/mnt/Koa/KoaDatasets_12152025/AK.SURV.csv")
d$died <- as.integer(d$Status.1 != "live")
d <- d[with(d, DBH.0>0 & Status.0=="live" & YIP>0 & is.finite(BYI) & is.finite(SDI.0)),]
d$exp <- pmax(d$YIP,0.1); d$RD <- d$SDI.0/500

cat("=== BYI distribution by origin (can plantations even span BYI?) ===\n")
cat("NATURAL  BYI:", paste(round(quantile(d$BYI[d$Planted==0],c(0,.1,.5,.9,1))),collapse=" / "),
    " n=",sum(d$Planted==0),"\n")
cat("PLANTED  BYI:", paste(round(quantile(d$BYI[d$Planted==1],c(0,.1,.5,.9,1))),collapse=" / "),
    " n=",sum(d$Planted==1),"\n")
cat("plantation records with BYI>399 (high tercile):", sum(d$Planted==1 & d$BYI>399),"\n")
cat("plantation deaths total:", sum(d$died[d$Planted==1]), " natural deaths:", sum(d$died[d$Planted==0]),"\n")

cat("\n=== annual mortality by origin x BYI tercile (ALL densities) ===\n")
d$bt <- cut(d$BYI, quantile(d$BYI,c(0,.33,.67,1),na.rm=T), include.lowest=T)
tab <- tapply(d$died,list(origin=d$Planted,byi=d$bt),sum)/tapply(d$exp,list(origin=d$Planted,byi=d$bt),sum)
print(round(tab,4))
cat("(cells with NA/0 = no plantation records there)\n")
cat("\ncounts:\n"); print(tapply(d$died,list(origin=d$Planted,byi=d$bt),length))

cat("\n=== Formal interaction test (Poisson, offset log exposure) ===\n")
m_main <- glm(died~Planted+I(log(BYI))+sqrt(SDI.0), poisson, offset=log(exp), data=d)
m_int  <- glm(died~Planted*I(log(BYI))+sqrt(SDI.0), poisson, offset=log(exp), data=d)
cat(sprintf(" main-effects AIC %.1f | with origin:BYI AIC %.1f (dAIC %.1f)\n",
    AIC(m_main),AIC(m_int),AIC(m_main)-AIC(m_int)))
cat(" interaction coef:", tryCatch(round(coef(m_int)["Planted:I(log(BYI))"],3),error=function(e)NA),
    " p=", tryCatch(signif(summary(m_int)$coef["Planted:I(log(BYI))","Pr(>|z|)"],3),error=function(e)NA),"\n")
cat(" (if SE is huge or coef extreme -> not identifiable due to confounding)\n")
print(round(summary(m_int)$coef,3))
