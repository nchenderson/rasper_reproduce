library(rasper)

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/RankDataIntegration/rasper_reproduce")
source("ExternalFunctions.R")

ProstateICI <- read.csv("DataAnalysis/MSKProstateICI.csv")
ProstateICI$ypseudo_3yr_cent <- ProstateICI$ypseudo_3yr - mean(ProstateICI$ypseudo_3yr)
Yct <- ProstateICI$ypseudo_3yr_cent

################################################
## Conventional Covariates
############################################

## 1. Duration of first-line chemotherapy
## 2. PSA prior to ICI (> 100 vs. <= 100)
## 3. ECOG performance status (< 1 vs. >= 1)
## 4. Presence of Visceral metastasis

#############################################
## Novel Covariates
############################################
## 1. MSI status
## 2. Presence of an alteration in a tumor-suppressor genes (TSG): TP53, PTEN, or RB1
## 3. Presence of a CDK12 alteration
## 4. tumor mutational burden.
risk_mod_form <- ypseudo_3yr_cent ~ psa + VISCERAL + ecog + chemo_dur + MSI_STATUS + TMB + mmrd + TSG + CDK12

lm_test <- lm(risk_mod_form, data=ProstateICI)

XX <- model.matrix(risk_mod_form, data=ProstateICI)
Xs <- scale(XX, scale=FALSE)[,-1]

## Define target vector of regression coefficients
beta.orig <- -c(log(2.09), log(1.63), log(1.91), log(1.57))
beta.tl <- c(beta.orig, rep(0, 6))
#beta.tl <- beta.orig
mu.ext <- as.numeric(Xs%*%beta.tl)

## The matrix RegCoefs will store the regression coefficients estimated by each model
nmethods <- 7
RegCoefs <- matrix(NA, nrow=ncol(Xs), ncol=nmethods)

## 1. OLS
ols_coefs <- as.numeric(solve(crossprod(Xs, Xs), crossprod(Xs, Yct)))
ols_fitted <- as.numeric(Xs%*%ols_coefs)
RegCoefs[,1] <- ols_coefs

## 2. Basic Ridge
ngrid <- 100
aseq <- c(0, exp(log(50) - log(100) + 0:ngrid*(log(100))/ngrid)) ## 200 leads to sum(coef^2) <= 0.2
loo_ridge <- LooCVNoninformativeRidge(y=Yct, X=Xs, alpha.seq=aseq)
alpha.ridge.best <- aseq[which.min(loo_ridge)]

direct_ridge <- NonInformativeRidge(y=Yct, X=Xs, alpha=alpha.ridge.best)
RegCoefs[,2] <- direct_ridge$coef

ridge_fitted <- direct_ridge$fitted.values

ngrid <- 20
## 3. Ranking Shrink (Kendall)
amax <- 50

aseq <- c(0, exp(log(amax) - log(100) + 0:ngrid*(log(100))/ngrid))

lmax <- exp(18)
nlamm <- 25
nlam <- nlamm + 1
bb <- 1000
lvals <- exp(log(lmax) - log(bb) + 0:nlamm*(log(bb))/nlamm)
CorVals <- CorValsS <- rep(0, nlam)
LamVals <- rep(0, nlam)
for(j in 1:nlam) {
  lseq <- lvals[j]
  tst_loo_kendall <- loocv.rasper(y=Yct, X=Xs, external.scores=mu.ext, 
                                  discrepancy="kendall", lambda.seq = lseq, 
                                  alpha.seq=aseq, maxiter=20)
  
  ind.best <- which.min(tst_loo_kendall$LooScore)
  alpha.best.kend <- tst_loo_kendall$alpha[ind.best]
  
  tst_loo_spear <- loocv.rasper(y=Yct, X=Xs, external.scores=mu.ext, 
                                discrepancy="spearman", lambda.seq = lseq, 
                                alpha.seq=aseq, maxiter=20)
  
  ind.best <- which.min(tst_loo_spear$LooScore)
  alpha.best.spear <- tst_loo_spear$alpha[ind.best]
  
  rshrink.kendall <- rasper(y=Yct, X=Xs, external.scores=mu.ext,
                            lambda = lseq, alpha = alpha.best.kend,
                            discrepancy="kendall", maxiter=500)
  
  rshrink.spear <- rasper(y=Yct, X=Xs, external.scores=mu.ext,
                          lambda = lseq, alpha = alpha.best.spear,
                          discrepancy="spearman", maxiter=500)
  
  rasp_kend_fitted <- rshrink.kendall$fitted.values
  rasp_spear_fitted <- rshrink.spear$fitted.values
  
  print(j)
  CorVals[j] <- cor(rasp_kend_fitted, mu.ext, method="kendall")
  CorValsS[j] <- cor(rasp_spear_fitted, mu.ext, method="kendall")
  LamVals[j] <- lseq
}


pdf("TablesFigures/Figure2.pdf",width=5.5, height=4.0)
par(mar=c(4.1, 4.1, 1.5, 0.5))
plot(log(LamVals), CorVals, xlab=expression(ln(lambda)), ylab=expression(paste("Kendall's ", tau)),
     main="Rank Association between the External Risk Model and RASPER", type="n",
     las=1, ylim=c(0.4, 1.0), cex.lab=1.2, cex.main=0.9)
points(log(LamVals), CorVals, pch=16)
points(log(LamVals), CorValsS, pch=16, col="red")
lines(log(LamVals), CorVals)
lines(log(LamVals), CorValsS, col="red")
legend("topleft", legend=c("RASPER (Spearman)", expression(paste("RASPER (Kendall's ", tau, ")"))), lwd=3,
       col=c("red", "black"), bty="n")
dev.off()
