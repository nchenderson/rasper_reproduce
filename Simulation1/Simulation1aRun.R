
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/RankDataIntegration/rankingshrink_reproduce")
source("ExternalFunctions.R")

library(rasper)
p <- 5
beta.true <- c(2, -1, 0.5, 0, 0)
######################################################
##. The only thing you need to change is nint if you want a different size
##. of the internal dataset
##. We considered nint=25, nint=50, and nint=100 so far
##################################################
nint <- 100

nreps <- 100 #(used for paper)

#nreps <- 3

BetaExt <- rbind(beta.true,
                 c(0.9, -0.1, 0.5, 0.25, -0.2),
                 c(3.6, -1.0, 2.5, 1.0, -0.8),
                 c(0.6, 0, -0.15, 0.3, 0.0),
                 c(3.5, -1.0, -1.3, -2.0, 1.0),
                 c(0.6, 0.1, -0.35, 0.2, -0.1),
                 c(2.8, -0.5, -1.5, -2.5, 1.5),
                 c(0.3, 0, -0.38, 0.35, -0.25),
                 c(1.5, 0.65, -0.5, 2.25, -1.5),
                 c(-0.5, -1, 0.3, 0, -0.1),
                 3*c(-0.5, -1, 0.3, 0, -0.1))

ngrid <- 20 ## number of grid points for cross-validation
nsettings <- nrow(BetaExt)
nmethods <- 8
TrueRankCorr <- MuDiscr <- matrix(NA, nreps, nsettings)
MSELM <- MSERPEN <- MSERPENAIC <- MSERPENK <- matrix(NA, nreps, nsettings)
MSEAug <- MSERidge <- MSEDTL <- MSEAngle <- matrix(NA, nreps, nsettings)
Results <- matrix(NA, nsettings, nmethods + 3)
CorMat <- matrix(NA, nsettings, 2)
TmpCorMat <- matrix(NA, nreps, 2)

set.seed(13572)

for(j in 1:nsettings) {
  beta.ext <- BetaExt[j,]
  for(h in 1:nreps) {
    X <- matrix(rnorm(nint*p), nrow=nint, ncol=p)
    Xnew <- matrix(rnorm(nint*p), nrow=nint, ncol=p)
    X <- scale(X)
    Xnew <- scale(Xnew)
    mu.true <- as.numeric(X%*%beta.true)
    mu.true.new <- as.numeric(Xnew%*%beta.true)
    Y <- X%*%beta.true + rnorm(nint, sd=5)
    Yprime <- X%*%beta.true + rnorm(nint, sd=5)

    mu.ext <- as.numeric(X%*%beta.ext)
    mu.ext.new <- as.numeric(X%*%beta.ext)
    external.ranking <- rank(as.numeric(mu.ext))/length(as.numeric(mu.ext))
    external.ranking.new <- rank(as.numeric(mu.ext.new))/length(as.numeric(mu.ext.new))
    Xr <- cbind(X, external.ranking)
    Xr.new <- cbind(Xnew, external.ranking.new)

    TrueRankCorr[h,j] <- cor(external.ranking, mu.true, method="spearman")
    MuDiscr[h,j] <- mean((mu.true - mu.ext)^2)

    ##################################################
    ### Run Different Methods
    ##################################################

    ## 1. OLS
    ols <- lm(Y ~ X - 1)
    ols.fitted <- ols$fitted.values
    ols.fitted <- as.numeric(Xnew%*%ols$coefficients)

    ## 2. Basic Ridge
    aseq <- c(0, exp(log(200) - log(100) + 0:ngrid*(log(100))/ngrid)) ## 200 leads to sum(coef^2) <= 0.2
    loo_ridge <- LooCVNoninformativeRidge(y=Y, X=X, alpha.seq=aseq)
    alpha.ridge.best <- aseq[which.min(loo_ridge)]

    direct_ridge <- NonInformativeRidge(y=Y, X=X, alpha=alpha.ridge.best)
    ridge.fitted <- direct_ridge$fitted.values
    ridge.fitted <- as.numeric(Xnew%*%direct_ridge$coef)

    ## 3. Distance TL
    aseq <- c(0, exp(log(400) - log(100) + 0:ngrid*(log(100))/ngrid)) ## 400 leads to sum(coef^2) <= 0.2
    loo_dtl <- LooCVInformativeRidge(y=Y, X=X, alpha.seq=aseq, beta.prior=beta.ext)
    alpha.dtl.best <- aseq[which.min(loo_dtl)]

    direct_dtl <- InformativeRidge(y=Y, X=X, alpha=alpha.dtl.best, beta.prior=beta.ext)
    dtl.fitted <- direct_dtl$fitted.values
    dtl.fitted <- as.numeric(Xnew%*%direct_dtl$coef)

    ## 3. Angle TL
    aseq <- c(0, exp(log(200) - log(100) + 0:ngrid*(log(100))/ngrid)) ## 200 leads to sum(coef^2) <= 0.2
    lseq <- c(0, exp(log(1000) - log(100) + 0:ngrid*(log(100))/ngrid)) ## 1000 leads to cor(muhat, mu.ext) >= 0.99 in most cases

    loo_angle <- LooCVAngleRidge(y=Y, X=X, alpha.seq=aseq, lambda.seq=lseq, beta.prior=beta.ext)
    alpha.angle.best <- loo_angle$alpha[which.min(loo_angle$LooScore)]
    lambda.angle.best <- loo_angle$lambda[which.min(loo_angle$LooScore)]

    direct_angle <- AngleRidge(y=Y, X=X, alpha=alpha.angle.best,
                               lambda=lambda.angle.best, beta.prior=beta.ext)
    angle.fitted <- direct_angle$fitted.values
    angle.fitted <- as.numeric(Xnew%*%direct_angle$coef)

    ## 4. RASPER (Spearman) with selection by AIC
    lmax <- 1000 ## lmax = 5000 gives rankcorr >= 0.94
    amax <- 500   ## amax = 500 gives sum(coef^2) <= 0.2
    lseq <- c(0, exp(log(lmax) - log(100) + 0:ngrid*(log(100))/ngrid))
    aseq <- c(0, exp(log(amax) - log(100) + 0:ngrid*(log(100))/ngrid))
    tst_loo <- loocv.rasper(y=Y, X=X, external.scores=mu.ext, lambda.seq = lseq,
                            alpha.seq=aseq, maxiter=20)

    aic.best <- which.min(tst_loo$aic)
    ind.best <- which.min(tst_loo$LooScore)
    alpha.best.rshrink <- tst_loo$alpha[ind.best]
    lambda.best.rshrink <- tst_loo$lambda[ind.best]
    alpha.best.aic <- tst_loo$alpha[aic.best]
    lambda.best.aic <- tst_loo$lambda[aic.best]

    rshrink <- rasper(y=Y, X=X, external.scores=mu.ext,
                      lambda = lambda.best.rshrink, alpha = alpha.best.rshrink,
                      discrepancy="spearman", maxiter=500)
    rshrink.fitted <- rshrink$fitted.values
    rshrink.fitted <- as.numeric(Xnew%*%rshrink$coef)

    rshrink.aic <- rasper(y=Y, X=X, external.scores=mu.ext,
                          lambda = lambda.best.aic, alpha = alpha.best.aic,
                          discrepancy="spearman", maxiter=500)
    rshrink.aic.fitted <- rshrink.aic$fitted.values
    rshrink.aic.fitted <- as.numeric(Xnew%*%rshrink.aic$coef)

    ## 5. RASPER (Kendall)
    lmax <- 500 ## lmax = 400 gives rankcorr >= 0.94
    amax <- 500   ## amax = 500 gives sum(coef^2) <= 0.2
    lseq <- c(0, exp(log(lmax) - log(100) + 0:ngrid*(log(100))/ngrid))
    aseq <- c(0, exp(log(amax) - log(100) + 0:ngrid*(log(100))/ngrid))
    tst_loo_kendall <- loocv.rasper(y=Y, X=X, external.scores=mu.ext, discrepancy="kendall",
                                    lambda.seq = lseq, alpha.seq=aseq, maxiter=20)

    kend.best <- which.min(tst_loo_kendall$LooScore)
    alpha.best.kendall <- tst_loo_kendall$alpha[kend.best]
    lambda.best.kendall <- tst_loo_kendall$lambda[kend.best]

    print(c(lambda.best.rshrink, lambda.best.kendall, alpha.best.rshrink, alpha.best.kendall))
    rshrink.kendall <- rasper(y=Y, X=X, external.scores=mu.ext,
                              lambda = lambda.best.kendall, alpha = alpha.best.kendall,
                              discrepancy="kendall", maxiter=500)
    rshrink.kendall.fitted <- rshrink.kendall$fitted.values
    rshrink.kendall.fitted <- as.numeric(Xnew%*%rshrink.kendall$coef)

    ## 6. Ranking Augmented Approach
    augm <- lm(Y ~ Xr - 1)
    augm.fitted <- augm$fitted.values
    augm.fitted <- as.numeric(Xr.new%*%augm$coefficients)

    #######################################
    ### Measure Performance
    #######################################

    MSELM[h,j] <- mean((ols.fitted - mu.true.new)^2)
    MSERidge[h,j] <- mean((ridge.fitted - mu.true.new)^2)
    MSEDTL[h,j] <- mean((dtl.fitted - mu.true.new)^2)
    MSEAngle[h,j] <- mean((angle.fitted - mu.true.new)^2)
    MSERPEN[h,j] <- mean((rshrink.fitted - mu.true.new)^2)
    MSEAug[h,j] <- mean((augm.fitted - mu.true.new)^2)
    MSERPENAIC[h,j] <- mean((rshrink.aic.fitted - mu.true.new)^2)
    MSERPENK[h,j] <- mean((rshrink.kendall.fitted - mu.true.new)^2)

    ## Use mu.true.new instead?

    corr.ridge<- cor(ridge.fitted, mu.ext, method="spearman")
    corr.rank <- cor(rshrink.fitted, mu.ext, method="spearman")
    TmpCorMat[h,] <- c(corr.ridge, corr.rank)
    print(c(h,j))
  }
  CorMat[j,] <- colMeans(TmpCorMat)
}
colnames(Results) <- c("Setting", "RankCorr", "MuDiscr", "Rank Pen", "Least Squares", "Ridge", "DTL", "ATL", "Augm",
                       "Rank Pen AIC", "Rank Pen Kendall")
Results[,1] <- 1:nsettings
Results[,2] <- round(colMeans(TrueRankCorr), 2)
Results[,3] <- round(colMeans(MuDiscr), 2)
Results[,4] <- colMeans(MSERPEN)
Results[,5] <- colMeans(MSELM)
Results[,6] <- colMeans(MSERidge)
Results[,7] <- colMeans(MSEDTL)
Results[,8] <- colMeans(MSEAngle)
Results[,9] <- colMeans(MSEAug)
Results[,10] <- colMeans(MSERPENAIC)
Results[,11] <- colMeans(MSERPENK)

if(nint==25) {
  fname <- "SimulationResults/Simulation1a_n25_tmp.RData"
} else if(nint==50) {
  fname <- "SimulationResults/Simulation1a_n50_tmp.RData"
} else if(nint==100) {
  fname <- "SimulationResults/Simulation1a_n100_tmp.RData"
}
save(Results, TrueRankCorr, MuDiscr, MSERPEN, MSELM, MSERidge, MSEDTL,
     MSEAngle, MSEAug, MSERPENAIC, MSERPENK, file=fname)
