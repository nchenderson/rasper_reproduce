setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/RankDataIntegration/rankingshrink_reproduce")
#setwd("~/Documents/RankDataIntegration/rankingshrink_reproduce")

library(rasper)
source("ExternalFunctions.R")


rtmix <- function(n, pp, df) {
  xx <- rep(NA, n)
  zz <- sample(0:1, size=n, replace=TRUE, prob=pp)
  xx[zz==0] <- rnorm(sum(zz==0))
  xx[zz==1] <- abs(rt(sum(zz==1), df=df))
  return(xx)
}

ff1 <- function(x) {
  ans <- plogis(x) - plogis(x - 1)
  return(ans)
}

ff2 <- function(x) {
  y <- exp(x)
  ans <- (0.5*(y - 2)^2)*(y < 7) + 12.5*(y > 7)
  return(ans)
}

ExternalRiskFunction <- function(X, theta.vec) {
  ans1 <- theta.vec[1]*ff1(X[,1]) + theta.vec[2]*ff2(X[,2]) + theta.vec[3]*ff1(X[,1])*ff2(X[,2])
  xthresh1 <- X[,3]*(X[,3] < 2)
  xthresh2 <- X[,3]*(X[,3] > 2)
  ans2 <- exp(-theta.vec[4]*xthresh1 + 10*theta.vec[4]*xthresh2 - theta.vec[5]*X[,5])
  ans <- (1 + ans1)*ans2
  return(ans)
}

nint <- 100
p <- 4
q <- 2
nXdraw <- 5

Btrans <- cbind(c(0.4, 0, 0.2, 0),
                c(0.25, 0, 0.6, 0.25))

ngrid <- 20 ## number of grid points for cross-validation
nsettings <- 10
nreps <- 100
nmethods <- 8
TrueRankCorr <- MuDiscr <- matrix(NA, nreps, nsettings)
MSELM <- MSERPEN <- MSERidge <- MSERPENM <- MSERPENK <- matrix(NA, nreps, nsettings)
MSEAug <- MSEDTL <- MSEAngle <- matrix(NA, nreps, nsettings)
MAELM <- MAERPEN <- MAERidge <- MAERPENM <- matrix(NA, nreps, nsettings)
MAEStack <- MAEDTL <- MAEAngle <- matrix(NA, nreps, nsettings)
Results <- matrix(NA, nsettings, nmethods + 3)

set.seed(13572)
ZBCov <- matrix(0, nrow=2, ncol=4)

BetaVals <- rbind(c(-0.6, 0.8, 1.5, -1.0, -1.2,0.0),
                  c(-0.25, 0.4, 0.7, -0.5, -0.5,0.0),
                  c(0.75, 1.6, 5.0, 0.0,1.5,0.0),
                  c(0.25, 0.4, 1.5, 0.0,0.5,0.0),
                  c(0.25, 0.0, 1.5, 0.0,0.0,0.0),
                  c(0.5, 0.0, 2.5, 0.0,0.0,0.0),
                  c(0.0, -5.0, 0.0, 0.0,0.0,0.0),
                  c(1.50, -2.5, 0.0, 0.0,0.0,0.0),
                  c(3.0, 0.0, 0.0, 0.0,0.0,0.0),
                  c(3.0, 0.0, 0.0, 0.0,0.0,0.0))
ThetaVals <- rbind(c(1.0, 0.5, 0.25, -0.2, 0.0),
                   c(1.0, 0.5, 0.25, -0.2, 0.0),
                   c(1.0, 0.5, 0.0, -0.2, 0.0),
                   c(1.0, 0.5, 0.0, -0.2, 0.0),
                   c(0.0, 0.0, 0.0, -0.2, 0.0),
                   c(0.0, 0.0, 0.0, -0.2, 0.0),
                   c(0.0, 1.0, 0.0, 0.0, 0.0),
                   c(1.0, 0.5, 0.0, 0.0, 0.0),
                   c(1.0, 0.0, 0.0, 0.0, 0.0),
                   c(6.0, 0.0, 0.0, 0.0, 0.0))

Z <- Znew <- matrix(NA, nrow=nint, ncol=p)
for(j in 1:nsettings) {
  theta.vecs <- ThetaVals[j,]
  beta.true <- BetaVals[j,]
  for(h in 1:nreps) {
    for(u in 1:p) {
      Z[,u] <- rtmix(nint, pp=c(0.75, 0.25), df=2)
      Znew[,u] <- rtmix(nint, pp=c(0.75, 0.25), df=2)
    }
    Z <- scale(Z)
    Znew <- scale(Znew)

    EpsB <- matrix(rnorm(nint*q), nrow=nint, ncol=q)
    EpsBNew <- matrix(rnorm(nint*q), nrow=nint, ncol=q)
    B <- Z%*%Btrans + EpsB
    Bnew <- Znew%*%Btrans + EpsBNew
    B <- scale(B)
    Bnew <- scale(Bnew)
    X <- cbind(Z, B)
    Xnew <- X

    ## sample covariance matrix between b_i and z_i
    ZBcov <- crossprod(B, Z)/nint
    mu.draw.marg <- Z%*%t(ZBcov)
    ## Get Cov matrix for drawing B in marginalized ranking pars.
    ZBZBt <- tcrossprod(ZBcov, ZBcov)
    VBdraw <- diag(c(1,1)) - ZBZBt
    LVBdraw <- try(t(chol(VBdraw)), silent=TRUE)
    if(class(LVBdraw)[1]=="try-error") {
      ZBcov <- 0.5*ZBcov
      ZBZBt <- tcrossprod(ZBcov, ZBcov)
      VBdraw <- diag(c(1,1)) - ZBZBt
      LVBdraw <- try(t(chol(VBdraw)), silent=TRUE)
      if(class(LVBdraw)[1]=="try-error") {
        LVBdraw <- diag(c(1,1))
      }
    }

    mu.true <- as.numeric(X%*%beta.true)
    mu.true.new <- as.numeric(Xnew%*%beta.true)
    mu.true <- mu.true - mean(mu.true)
    mu.true.new <- mu.true.new - mean(mu.true)
    Y <- mu.true + rnorm(nint, sd=5)

    mu.ext <- ExternalRiskFunction(X, theta.vecs)
    mu.ext.new <- ExternalRiskFunction(Xnew, theta.vecs)
    external.ranking <- rank(as.numeric(mu.ext))/length(as.numeric(mu.ext))
    external.ranking.new <- rank(as.numeric(mu.ext.new))/length(as.numeric(mu.ext.new))
    Xr <- cbind(X, external.ranking)
    Xr.new <- cbind(Xnew, external.ranking.new)

    TrueRankCorr[h,j] <- cor(external.ranking, mu.true, method="spearman")
    MuDiscr[h,j] <- mean((mu.true - mu.ext)^2)

    ##################################################
    ### Run Different Methods
    ##################################################

    ## First, define target vector for DTL and ATL:
    beta.zhat <- lm(mu.ext ~ X)$coef[-1]
    beta.tl <- beta.zhat
    print(round(beta.tl, 3))

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
    loo_dtl <- LooCVInformativeRidge(y=Y, X=X, alpha.seq=aseq, beta.prior=beta.tl)
    alpha.dtl.best <- aseq[which.min(loo_dtl)]

    direct_dtl <- InformativeRidge(y=Y, X=X, alpha=alpha.dtl.best, beta.prior=beta.tl)
    dtl.fitted <- direct_dtl$fitted.values
    dtl.fitted <- as.numeric(Xnew%*%direct_dtl$coef)

    ## 3. Angle TL
    aseq <- c(0, exp(log(200) - log(100) + 0:ngrid*(log(100))/ngrid)) ## 200 leads to sum(coef^2) <= 0.2
    lseq <- c(0, exp(log(1000) - log(100) + 0:ngrid*(log(100))/ngrid)) ## 1000 leads to cor(muhat, mu.ext) >= 0.99 in most cases

    loo_angle <- LooCVAngleRidge(y=Y, X=X, alpha.seq=aseq, lambda.seq=lseq, beta.prior=beta.tl)
    alpha.angle.best <- loo_angle$alpha[which.min(loo_angle$LooScore)]
    lambda.angle.best <- loo_angle$lambda[which.min(loo_angle$LooScore)]

    direct_angle <- AngleRidge(y=Y, X=X, alpha=alpha.angle.best,
                               lambda=lambda.angle.best, beta.prior=beta.tl)
    angle.fitted <- direct_angle$fitted.values
    angle.fitted <- as.numeric(Xnew%*%direct_angle$coef)

    ## 4. RASPER Spearman
    lmax <- 1000 ## lmax = 5000 gives rankcorr >= 0.94
    amax <- 500   ## amax = 500 gives sum(coef^2) <= 0.2
    lseq <- c(0, exp(log(lmax) - log(100) + 0:ngrid*(log(100))/ngrid))
    aseq <- c(0, exp(log(amax) - log(100) + 0:ngrid*(log(100))/ngrid))
    tst_loo <- loocv.rasper(y=Y, X=X, external.scores=mu.ext,
                            lambda.seq = lseq, alpha.seq=aseq, maxiter=20)

    aic.best <- which.min(tst_loo$aic)
    ind.best <- which.min(tst_loo$LooScore)
    alpha.best.rshrink <- tst_loo$alpha[ind.best]
    lambda.best.rshrink <- tst_loo$lambda[ind.best]
    alpha.best.aic <- tst_loo$alpha[aic.best]
    lambda.best.aic <- tst_loo$lambda[aic.best]


    rshrink <- rasper(y=Y, X=X, external.scores=mu.ext, lambda = lambda.best.rshrink,
                      alpha = alpha.best.rshrink, discrepancy="spearman",
                      maxiter=500)
    rshrink.fitted <- rshrink$fitted.values
    rshrink.fitted <- as.numeric(Xnew%*%rshrink$coef)

    ## 5. RASPER (Kendall)
    lmax <- 700 ## lmax = 700 gives about maximum rank correlatoin
    amax <- 500   ## amax = 500 gives sum(coef^2) <= 0.2
    lseq <- c(0, exp(log(lmax) - log(100) + 0:ngrid*(log(100))/ngrid))
    aseq <- c(0, exp(log(amax) - log(100) + 0:ngrid*(log(100))/ngrid))
    tst_loo_kendall <- loocv.rasper(y=Y, X=X, external.scores=mu.ext, discrepancy="kendall",
                                    lambda.seq = lseq, alpha.seq=aseq, maxiter=20)

    kend.best <- which.min(tst_loo_kendall$LooScore)
    alpha.best.kendall <- tst_loo_kendall$alpha[kend.best]
    lambda.best.kendall <- tst_loo_kendall$lambda[kend.best]

    #print(c(lambda.best.rshrink, lambda.best.kendall, alpha.best.rshrink, alpha.best.kendall))
    rshrink.kendall <- rasper(y=Y, X=X, external.scores=mu.ext,
                              lambda = lambda.best.kendall, alpha = alpha.best.kendall,
                              discrepancy="kendall", maxiter=500)
    rshrink.kendall.fitted <- rshrink.kendall$fitted.values
    rshrink.kendall.fitted <- as.numeric(Xnew%*%rshrink.kendall$coef)

    ## 6. RASPER with marginalized parameters

    Xlist <- list()
    Xlist[[1]] <- X
    for(u in 1:nXdraw) {
      eps.b.marg <- matrix(rnorm(nint*q), nrow=q, ncol=nint)
      Btmpp <- mu.draw.marg + t(LVBdraw%*%eps.b.marg)
      Xlist[[u+1]] <- cbind(Z, Btmpp)
    }

    rshrink.marg <- rasper(y=Y, X=X, Xlist=Xlist, external.scores=mu.ext,
                           lambda = lambda.best.rshrink, alpha = alpha.best.rshrink,
                           discrepancy="spearman", maxiter=500)
    rshrink.marg.fitted <- rshrink.marg$fitted.values
    rshrink.marg.fitted <- as.numeric(Xnew%*%rshrink.marg$coef)


    ## 5. Ranking Augmentation
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
    MSERPENM[h,j] <- mean((rshrink.marg.fitted - mu.true.new)^2)
    MSERPENK[h,j] <- mean((rshrink.kendall.fitted - mu.true.new)^2)

    MAELM[h,j] <- median(abs(ols.fitted - mu.true.new))
    MAERidge[h,j] <- median(abs(ridge.fitted - mu.true.new))
    MAEDTL[h,j] <- median(abs(dtl.fitted - mu.true.new))
    MAEAngle[h,j] <- median(abs(angle.fitted - mu.true.new))
    MAERPEN[h,j] <- median(abs(rshrink.fitted - mu.true.new))
    #MAEStack[h,j] <- median(abs(stacked.fitted - mu.true.new))
    MAERPENM[h,j] <- median(abs(rshrink.marg.fitted - mu.true.new))

    ## Use mu.true.new instead?
    print(c(h,j))
  }
}

colnames(Results) <- c("Setting", "RankCorr", "MuDiscr", "Rank Pen", "Least Squares",
                       "Ridge", "DTL", "ATL", "Stacked",
                       "Rank Pen Marg","Rank Pen Kendall")

Results[,1] <- 1:nsettings
Results[,2] <- round(colMeans(TrueRankCorr), 2)
Results[,3] <- round(colMeans(MuDiscr), 2)
Results[,4] <- colMeans(MSERPEN)
Results[,5] <- colMeans(MSELM)
Results[,6] <- colMeans(MSERidge)
Results[,7] <- colMeans(MSEDTL)
Results[,8] <- colMeans(MSEAngle)
Results[,9] <- colMeans(MSEAug)
Results[,10] <- colMeans(MSERPENM)
Results[,11] <- colMeans(MSERPENK)

Results <- Results[order(-Results[,2]),]
Results[,1] <- 1:10
round(Results, 2)

if(nint==25) {
  fname <- "SimulationResults/Simulation2a_n25_tmp.RData"
} else if(nint==50) {
  fname <- "SimulationResults/Simulation2a_n50_tmp.RData"
} else if(nint==100) {
  fname <- "SimulationResults/Simulation2a_n100_tmp.RData"
}

save(Results, TrueRankCorr, MuDiscr, MSERPEN, MSELM, MSERidge, MSEDTL,
     MSEAngle, MSEAug, MSERPENK, MSERPENM, MAERPEN, MAELM, MAERidge, MAEDTL,
     MAEAngle, MAEStack, MAERPENM, file=fname)



