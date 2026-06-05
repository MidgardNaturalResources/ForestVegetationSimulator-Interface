# extract_ingrowth_psp.R -- derive koa ingrowth (new trees between measurements)
# from AK_PSP_data_compiled.csv. Produces results/koa_ingrowth_obs.csv used by the
# RD-only exploratory fits (fit_ingrowth.R / fit_ingrowth_rd.R). The canonical
# ingrowth model uses AK.HT.csv (fit_ingrowth_byi.R) which carries BYI/origin.
d <- read.csv("/sessions/elegant-lucid-cannon/mnt/Koa/AK_PSP_data_compiled.csv")
d <- d[d$DBH_cm>0,]
psp <- split(d, d$PSP); ing <- data.frame()
for(p in names(psp)){
  g <- psp[[p]]; seqs <- sort(unique(g$MeasSeq)); if(length(seqs)<2) next
  for(i in 2:length(seqs)){
    prev <- g$Tree_no[g$MeasSeq==seqs[i-1]]; cur <- g[g$MeasSeq==seqs[i],]
    new <- cur[!(cur$Tree_no %in% prev),]
    yrs <- (unique(cur$Yfloor)[1]-unique(g$Yfloor[g$MeasSeq==seqs[i-1]])[1])
    if(is.na(yrs)||yrs<=0) next
    tph_prev <- sum(g$TPA[g$MeasSeq==seqs[i-1]])
    ba_prev  <- sum((g$DBH_cm[g$MeasSeq==seqs[i-1]]^2*0.00007854)*g$TPA[g$MeasSeq==seqs[i-1]])
    ing <- rbind(ing, data.frame(PSP=p, seq=seqs[i], yrs=yrs, n_new=nrow(new),
      new_tph=sum(new$TPA), ann_ingrowth_tph=sum(new$TPA)/yrs,
      BAPH=ba_prev, TPH=tph_prev, minDBH_new=ifelse(nrow(new)>0,min(new$DBH_cm),NA)))
  }
}
write.csv(ing, "/sessions/elegant-lucid-cannon/mnt/outputs/koa_harness/results/koa_ingrowth_obs.csv", row.names=FALSE)
cat("Wrote koa_ingrowth_obs.csv:", nrow(ing), "plot-periods,", sum(ing$n_new>0), "with ingrowth\n")
