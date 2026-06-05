# fit_ingrowth_byi.R -- reconstruct koa ingrowth from AK.HT.csv (has BYI, SDI,
# Origin, pBA.AK) and test RD, BYI, and BYI x RD interaction. koa-only.
h <- read.csv("/sessions/elegant-lucid-cannon/mnt/Koa/KoaDatasets_12152025/AK.HT.csv")
h$key <- paste(h$Data,h$Install,h$Plot,sep="|")
h$tid <- paste(h$key,h$Tree,sep="|")
yrcol <- if("Year"%in%names(h)) "Year" else "Age"
ing <- list()
for(k in unique(h$key)){
  g <- h[h$key==k,]; ms <- sort(unique(g$Measure))
  if(length(ms)<2) next
  for(i in 2:length(ms)){
    prev <- g[g$Measure==ms[i-1],]; cur <- g[g$Measure==ms[i],]
    yrs <- (cur[[yrcol]][1]-prev[[yrcol]][1]); if(is.na(yrs)||yrs<=0||yrs>30) next
    new <- cur[!(cur$tid %in% prev$tid),]
    new <- new[!is.na(new$Status) & new$Status %in% c("live","Live","L",1) | TRUE,] # keep all new
    s <- prev[1,]
    ing[[length(ing)+1]] <- data.frame(key=k, planted=as.integer(s$Origin=="Planted"|s$Origin==1),
      yrs=yrs, n_new=nrow(new), ing_tph=sum(new$EXPF,na.rm=TRUE)/yrs,
      BAPH=s$BAPH, SDI=s$SDI, RD=s$SDI/500, BYI=s$BYI, pBA=s$pBA.AK)
  }
}
d <- do.call(rbind, ing); d <- d[is.finite(d$RD)&is.finite(d$BYI),]
cat("plot-periods:",nrow(d)," %occur:",round(mean(d$ing_tph>0),2),
    " mean ingrowth:",round(mean(d$ing_tph),1)," | planted n=",sum(d$planted),"\n")
cat("pBA.AK (pct koa BA): median",round(median(d$pBA,na.rm=T),2)," IQR",
    paste(round(quantile(d$pBA,c(.25,.75),na.rm=T),2),collapse="-"),"\n")
cat("BYI range:",paste(round(range(d$BYI),0),collapse="-")," by origin: nat",
    round(median(d$BYI[d$planted==0]),0)," plt",round(median(d$BYI[d$planted==1]),0),"\n\n")

cat("Mean ingrowth by RD x BYI tercile:\n")
d$rb<-cut(d$RD,c(0,.3,.6,3),labels=c("RD<.3","RD.3-.6","RD>.6"))
d$bb<-cut(d$BYI,quantile(d$BYI,c(0,.33,.67,1),na.rm=T),include.lowest=T,labels=c("BYIlo","BYImid","BYIhi"))
print(round(tapply(d$ing_tph,list(d$rb,d$bb),mean,na.rm=TRUE),1))

cat("\n=== Models (quasipoisson, log link) ===\n")
m_rd  <- glm(ing_tph~RD, quasipoisson, data=d)
m_rdb <- glm(ing_tph~RD+I(log(BYI)), quasipoisson, data=d)
m_int <- glm(ing_tph~RD*I(log(BYI)), quasipoisson, data=d)
m_intO<- glm(ing_tph~RD*I(log(BYI))+planted, quasipoisson, data=d)
cat("RD only        coef:",paste(round(coef(m_rd),4),collapse=", "),"\n")
cat("RD+BYI         coef:",paste(round(coef(m_rdb),4),collapse=", "),
    " BYI p=",signif(summary(m_rdb)$coef["I(log(BYI))","Pr(>|t|)"],3),"\n")
cat("RD*BYI         coef:",paste(round(coef(m_int),4),collapse=", "),"\n")
cat("  interaction RD:log(BYI) p=",signif(summary(m_int)$coef["RD:I(log(BYI))","Pr(>|t|)"],3),"\n")
cat("RD*BYI+planted planted p=",tryCatch(signif(summary(m_intO)$coef["planted","Pr(>|t|)"],3),error=function(e)NA),"\n")

cat("\nBehavior: E[ingrowth] at BYI lo/mid/hi across RD (RD*BYI model)\n")
gr<-expand.grid(RD=c(.1,.3,.5,.7),BYI=quantile(d$BYI,c(.2,.5,.8)))
gr$E<-predict(m_int,newdata=gr,type="response")
gr$BYI<-round(gr$BYI); print(round(gr,1))
write.csv(d,"/sessions/elegant-lucid-cannon/mnt/outputs/koa_harness/results/koa_ingrowth_byi_obs.csv",row.names=FALSE)
co<-coef(m_int); con<-file("/sessions/elegant-lucid-cannon/mnt/outputs/koa_harness/results/ingrowth_byi_coef.csv","w")
cat("term,estimate\n",file=con); for(t in names(co)) cat(sprintf("%s,%.8f\n",t,co[t]),file=con); close(con)
cat("\nWrote koa_ingrowth_byi_obs.csv, ingrowth_byi_coef.csv\n")
