# fit_ingrowth_rd.R -- koa-only annualized ingrowth as a single function of
# relative density (RD = SDI/SDImax). Expected annual ingrowth (trees/ha/yr);
# the mean already absorbs the 74% zeros. Annualized by construction (rate/yr).
d <- read.csv("/sessions/elegant-lucid-cannon/mnt/outputs/koa_harness/results/koa_ingrowth_obs.csv")
d <- d[is.finite(d$BAPH) & is.finite(d$TPH) & d$yrs>0 & d$TPH>0,]
d$QMD <- sqrt(d$BAPH/(0.00007854*d$TPH))
d$SDI <- d$TPH*(d$QMD/25)^1.6
d$RD  <- d$SDI/500
d$rate <- d$ann_ingrowth_tph

cat("n=",nrow(d)," RD range",round(min(d$RD),2),"-",round(max(d$RD),2),
    " mean ingrowth",round(mean(d$rate),2),"trees/ha/yr\n\n")
cat("Observed mean annual ingrowth by RD band:\n")
d$rb <- cut(d$RD, c(0,.2,.4,.6,.85,5), labels=c("<.2",".2-.4",".4-.6",".6-.85",">.85"))
for(l in levels(d$rb)){s<-d[!is.na(d$rb)&d$rb==l,]
  cat(sprintf("  RD %-7s n=%3d  mean %5.1f  %%occur %2.0f\n",l,nrow(s),mean(s$rate),100*mean(s$rate>0)))}

# --- Form A: exponential decline  E = a*exp(-b*RD) ---
mA <- glm(rate ~ RD, family=quasipoisson(link="log"), data=d)
a<-exp(coef(mA)[1]); b<- -coef(mA)[2]
cat(sprintf("\nForm A (exp):  E = %.2f * exp(-%.2f * RD)\n", a, b))

# --- Form B: linear decline to zero at RDc  E = I0*max(0,1-RD/RDc) ---
# grid search RDc, I0 by least squares on observed mean
rss<-function(par){I0<-par[1];RDc<-par[2]; pred<-I0*pmax(0,1-d$RD/RDc); sum((d$rate-pred)^2)}
opt<-optim(c(20,0.9), rss, method="L-BFGS-B", lower=c(1,0.5), upper=c(80,2))
I0<-opt$par[1]; RDc<-opt$par[2]
cat(sprintf("Form B (linear): E = %.1f * max(0, 1 - RD/%.2f)\n", I0, RDc))

# --- compare predicted vs observed by RD band ---
cat("\nPredicted vs observed by RD band:\n  RD       obs   formA  formB\n")
for(l in levels(d$rb)){s<-d[!is.na(d$rb)&d$rb==l,]; rdm<-mean(s$RD)
  cat(sprintf("  %-7s %5.1f  %5.1f  %5.1f\n",l,mean(s$rate),
      a*exp(-b*rdm), I0*max(0,1-rdm/RDc)))}
cat(sprintf("\nRMSE  A %.2f  B %.2f\n",
  sqrt(mean((d$rate-a*exp(-b*d$RD))^2)),
  sqrt(mean((d$rate-I0*pmax(0,1-d$RD/RDc))^2))))

# write chosen coefficients (both, pick in projector)
con<-file("/sessions/elegant-lucid-cannon/mnt/outputs/koa_harness/results/ingrowth_rd_coef.csv","w")
cat("form,par,value\n",file=con)
cat(sprintf("expA,a,%.6f\nexpA,b,%.6f\n",a,b),file=con)
cat(sprintf("linB,I0,%.6f\nlinB,RDc,%.6f\n",I0,RDc),file=con)
close(con)
cat("\nWrote ingrowth_rd_coef.csv\n")
