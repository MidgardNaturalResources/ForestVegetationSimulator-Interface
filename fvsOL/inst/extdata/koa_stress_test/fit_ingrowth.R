# fit_ingrowth.R -- koa ingrowth equations following Li, Weiskittel & Kershaw
# (2011). Two-stage / zero-inflated structure: 74% of plot-periods have zero
# ingrowth. Test alternatives, pick by AIC + logical behavior.
d <- read.csv("/sessions/elegant-lucid-cannon/mnt/outputs/koa_harness/results/koa_ingrowth_obs.csv")
d <- d[is.finite(d$BAPH) & is.finite(d$TPH) & d$yrs>0,]
d$occur <- as.integer(d$n_new>0)
d$rate  <- d$ann_ingrowth_tph                      # trees/ha/yr
cat("n plot-periods:",nrow(d)," occurrence:",round(mean(d$occur),3),
    " mean rate:",round(mean(d$rate),2)," (nonzero mean:",round(mean(d$rate[d$occur==1]),1),")\n\n")

# --- Candidate 1: STATIC (constant) ---
static_mean <- mean(d$rate)
cat("C1 STATIC: E[ingrowth] =",round(static_mean,2),"trees/ha/yr (constant)\n\n")

# --- Candidate 2: Poisson on counts, offset log(yrs) [Li et al. eq.7 style] ---
p2 <- glm(n_new ~ BAPH + TPH, family=poisson, offset=log(yrs), data=d)
cat("C2 POISSON ln(lambda)=b0+b1*BAPH+b2*TPH (offset logYIP):\n")
print(round(coef(p2),5)); cat(" AIC",round(AIC(p2),1),"\n\n")

# --- Candidate 3: Negative binomial (overdispersion) ---
nb_ok <- requireNamespace("MASS", quietly=TRUE)
if(nb_ok){
  p3 <- MASS::glm.nb(n_new ~ BAPH + TPH + offset(log(yrs)), data=d)
  cat("C3 NEG BINOMIAL:\n"); print(round(coef(p3),5)); cat(" AIC",round(AIC(p3),1)," theta",round(p3$theta,3),"\n\n")
}

# --- Candidate 4: TWO-STAGE (Li et al. core): logistic occurrence x amount ---
occ <- glm(occur ~ BAPH + TPH, family=binomial, data=d)
amt <- glm(rate ~ BAPH + TPH, family=Gamma(link="log"), data=d[d$occur==1,])
cat("C4 TWO-STAGE:\n occurrence logit(P)=",paste(round(coef(occ),5),collapse=", "),"\n")
cat(" amount|occur log(E)=",paste(round(coef(amt),5),collapse=", "),"\n")
# expected = P(occur)*E[amount|occur]
d$pP <- predict(occ,type="response")
d$pA <- predict(amt,newdata=d,type="response")
d$pExp <- d$pP*d$pA
cat(" two-stage E[ingrowth]: mean",round(mean(d$pExp),2)," (obs mean",round(mean(d$rate),2),")\n\n")

# --- behavior: predicted ingrowth vs BAPH (should DECLINE; canopy closure) ---
cat("Logical-behavior check (TPH=400): expected annual ingrowth TPH vs BAPH\n")
grid <- data.frame(BAPH=c(2,5,10,20,30,45,60), TPH=400, yrs=1)
grid$C2 <- predict(p2,newdata=grid,type="response")
grid$pP <- predict(occ,newdata=grid,type="response"); grid$pA <- predict(amt,newdata=grid,type="response")
grid$C4 <- grid$pP*grid$pA
print(round(grid[,c("BAPH","C2","C4")],2))

# write chosen (two-stage) coefficients
co_occ<-coef(occ); co_amt<-coef(amt)
con<-file("/sessions/elegant-lucid-cannon/mnt/outputs/koa_harness/results/ingrowth_coef.csv","w")
cat("part,term,estimate\n",file=con)
for(t in names(co_occ)) cat(sprintf("occur,%s,%.8f\n",t,co_occ[t]),file=con)
for(t in names(co_amt)) cat(sprintf("amount,%s,%.8f\n",t,co_amt[t]),file=con)
close(con)
cat("\nWrote ingrowth_coef.csv (two-stage)\n")
