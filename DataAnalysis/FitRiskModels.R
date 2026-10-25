
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


#risk_mod_form <- ypseudo_3yr_cent ~ psa + VISCERAL + ecog + chemo_dur
risk_mod_form <- ypseudo_3yr_cent ~ psa + VISCERAL + ecog + chemo_dur + MSI_STATUS + TMB + mmrd + TSG + CDK12

lm_test <- lm(risk_mod_form, data=ProstateICI)

XX <- model.matrix(risk_mod_form, data=ProstateICI)
Xs <- scale(XX, scale=FALSE)[,-1]

## Define target vector of regression coefficients
beta.orig <- -c(log(2.09), log(1.63), log(1.91), log(1.57))
beta.tl <- c(beta.orig, rep(0, 6))
#beta.tl <- beta.orig
mu.ext <- as.numeric(Xs%*%beta.tl)
#cor(mu.ext, lm_test$fitted.values, method="spearman")

## The matrix RegCoefs will store the regression coefficients estimated by each model
nmethods <- 7
RegCoefs <- matrix(NA, nrow=ncol(Xs), ncol=nmethods)
FitValMat <- matrix(NA, nrow=nrow(Xs), ncol=nmethods)

## 1. OLS

ols_coefs <- as.numeric(solve(crossprod(Xs, Xs), crossprod(Xs, Yct)))
ols_fitted <- as.numeric(Xs%*%ols_coefs)
RegCoefs[,1] <- ols_coefs
FitValMat[,1] <- ols_fitted

## 2. Basic Ridge
ngrid <- 100
aseq <- c(0, exp(log(50) - log(100) + 0:ngrid*(log(100))/ngrid)) ## 200 leads to sum(coef^2) <= 0.2
loo_ridge <- LooCVNoninformativeRidge(y=Yct, X=Xs, alpha.seq=aseq)
alpha.ridge.best <- aseq[which.min(loo_ridge)]

direct_ridge <- NonInformativeRidge(y=Yct, X=Xs, alpha=alpha.ridge.best)
RegCoefs[,2] <- direct_ridge$coef
FitValMat[,2] <- direct_ridge$fitted.values

## 3. Distance TL
## Need to refine this target
aseq <- c(0, exp(log(50) - log(100) + 0:ngrid*(log(100))/ngrid)) ## 400 leads to sum(coef^2) <= 0.2
loo_dtl <- LooCVInformativeRidge(y=Yct, X=Xs, alpha.seq=aseq, beta.prior=beta.tl)

#plot(aseq, loo_dtl)

alpha.dtl.best <- aseq[which.min(loo_dtl)]
direct_dtl <- InformativeRidge(y=Yct, X=Xs, alpha=alpha.dtl.best, beta.prior=beta.tl)

RegCoefs[,3] <- direct_dtl$coef
FitValMat[,3] <- direct_dtl$fitted.values

## 3. Angle TL
lmax <- exp(12)
aseq <- c(0, exp(log(50) - log(100) + 0:ngrid*(log(100))/ngrid))
lseq <- c(0, exp(log(lmax) - log(100) + 0:ngrid*(log(100))/ngrid))

loo_angle <- LooCVAngleRidge(y=Yct, X=Xs, alpha.seq=aseq, lambda.seq=lseq, beta.prior=beta.tl)
alpha.angle.best <- loo_angle$alpha[which.min(loo_angle$LooScore)]
lambda.angle.best <- loo_angle$lambda[which.min(loo_angle$LooScore)]

direct_angle <- AngleRidge(y=Yct, X=Xs, alpha=alpha.angle.best,
                           lambda=lambda.angle.best, beta.prior=beta.tl)
RegCoefs[,4] <- direct_angle$coef
FitValMat[,4] <- direct_angle$fitted.values

## 5. Augmented Approach
X.ext <- rank(mu.ext)
Xr <- cbind(Xs, X.ext)
stacked <- lm(Yct ~ Xr - 1)
augm_fit <- as.numeric(Xr%*%stacked$coefficients)
RegCoefs[,5] <- stacked$coefficients[1:nrow(RegCoefs)]
FitValMat[,5] <- augm_fit

## 6. Ranking Shrink (Kendall)
ngrid <- 50
lmax <- exp(15)
amax <- 50
lseq <- c(0, exp(log(lmax) - log(100) + 0:ngrid*(log(100))/ngrid))
aseq <- c(0, exp(log(amax) - log(100) + 0:ngrid*(log(100))/ngrid))
tst_loo_kendall <- loocv.rasper(y=Yct, X=Xs, external.scores=mu.ext, 
                                discrepancy="kendall", lambda.seq = lseq, 
                                alpha.seq=aseq, maxiter=20)

kend.best <- which.min(tst_loo_kendall$LooScore)
alpha.best.kendall <- tst_loo_kendall$alpha[kend.best]
lambda.best.kendall <- tst_loo_kendall$lambda[kend.best]

rshrink.kendall <- rasper(y=Yct, X=Xs, external.scores=mu.ext,
                          lambda = lambda.best.kendall, alpha = alpha.best.kendall,
                          discrepancy="kendall", maxiter=500)

RegCoefs[,6] <- rshrink.kendall$coef
FitValMat[,6] <- rshrink.kendall$fitted.values

####

## 7. Ranking Shrink (Spearman)
lmax <- exp(15) ## lmax = 5000 gives rankcorr >= 0.94
amax <- 50   ## amax = 500 gives sum(coef^2) <= 0.2
lseq <- c(0, exp(log(lmax) - log(100) + 0:ngrid*(log(100))/ngrid))
aseq <- c(0, exp(log(amax) - log(100) + 0:ngrid*(log(100))/ngrid))

tst_loo_spear <- loocv.rasper(y=Yct, X=Xs, external.scores=mu.ext, 
                              lambda.seq = lseq, alpha.seq=aseq, 
                              maxiter=20)

spear.best <- which.min(tst_loo_spear$LooScore)
alpha.best.spear <- tst_loo_spear$alpha[spear.best]
lambda.best.spear <- tst_loo_spear$lambda[spear.best]

rshrink.spear <- rasper(y=Yct, X=Xs, external.scores=mu.ext,
                        lambda = lambda.best.spear, 
                        alpha = alpha.best.spear,
                        discrepancy="spearman", maxiter=500)

RegCoefs[,7] <- rshrink.spear$coef
FitValMat[,7] <- rshrink.spear$fitted.values

colnames(RegCoefs) <- c("OLS", "Ridge", "DTL", "ATL", "Stack", "RASPERK", "RASPERS")
rownames(RegCoefs) <- colnames(Xs)
RegCoefs <- RegCoefs[c(1, 2, 3, 4, 5, 6, 8, 7, 9, 10),]

## Save the regression coefficients and fitted values: 
save(RegCoefs, file="DataAnalysis/EstimatedRegressionCoefficients.RData")
save(FitValMat, file="DataAnalysis/AllFittedValues.RData")

