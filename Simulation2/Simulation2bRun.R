
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/RankDataIntegration/rankingshrink_reproduce")
source("ExternalFunctions.R")

library(rankingshrink)
p <- 4
q <- 2
beta.true <- c(2, -1, 0.5, 0, -0.2, 1.0)
nint <- 50

nreps <- 5000
nXdraw <- 0

ff1 <- function(x) {return(x/4)}
ff2 <- function(x) {return(0.5*(x > 0) - 0.5*(x <= 0)) }
ff3 <- function(x) {return(0.1*(x^3))}
ff4 <- function(x) {return((1/(1 + exp(-4*x))) + (- 0.1*x^2)*(x > 0))}
ff5 <- function(x) {return(-x/8)}

ExtRiskScoreFn <- function(Z, theta.vec) {
  Znew <- Z*theta.vec
  add.term <- ff1(Znew[,1]) + ff2(Znew[,2]) + ff3(Znew[,3]) + ff4(Znew[,4])
  int.term1 <- ff1(Znew[,1])*ff2(Znew[,2]) + ff1(Znew[,1])*ff3(Znew[,3])
  int.term2 <- ff2(Znew[,2])*ff3(Znew[,3]) + ff2(Znew[,2])*ff4(Znew[,4])
  int.term3 <- ff3(Znew[,3])*ff4(Znew[,4]) + ff1(Znew[,1])*ff4(Znew[,4])

  ans <- add.term + int.term1 + int.term2 + int.term3
  return(ans)
}

IntRiskScoreFn <- function(X, gamma.vec) {
  ## fill in
  Xnew <- X*gamma.vec
  add.term <- ff1(Xnew[,1]) + ff2(Xnew[,2]) + ff3(Xnew[,3]) + ff4(Xnew[,4]) + ff5(Xnew[,5])
}

Btrans <- cbind(c(0.4, 0, 0, 0),
                c(0.25, 0, 0.5, 0.2))

xx <- seq(-1, 1, length.out=8)
Beta.grid <- as.matrix(expand.grid(xx, xx, xx, xx, xx, xx, xx))

ngrid <- nrow(Beta.grid)
TrueRankCorr <- MuDiscr <- rep(NA, ngrid)
for(k in 1:ngrid) {
  beta.ext <- Beta.grid[k,]

  Z <- matrix(rnorm(nint*p), nrow=nint, ncol=p)
  Znew <- matrix(rnorm(nint*p), nrow=nint, ncol=p)
  Z <- scale(Z)
  Znew <- scale(Znew)
  EpsB <- matrix(rnorm(nint*q), nrow=nint, ncol=q)
  EpsBNew <- matrix(rnorm(nint*q), nrow=nint, ncol=q)

  B <- Z%*%Btrans + EpsB
  Bnew <- Znew%*%Btrans + EpsBNew
  B <- scale(B)
  Bnew <- scale(Bnew)
  X <- cbind(Z, B)
  Xnew <- cbind(Znew, Bnew)

  mu.true <- as.numeric(X%*%beta.true)
  mu.true.new <- as.numeric(Xnew%*%beta.true)
  Y <- X%*%beta.true + rnorm(nint, sd=5)
  Yprime <- X%*%beta.true + rnorm(nint, sd=5)

  mu.ext <- ExtRiskScoreFn(Z, beta.ext)
  external.ranking <- rank(as.numeric(mu.ext))/length(as.numeric(mu.ext))

  TrueRankCorr[k] <- cor(external.ranking, mu.true, method="spearman")
  MuDiscr[k] <- mean((mu.true - mu.ext)^2)
  if(k%%100000 == 0) {
    print(c(k, round(k/nrow(Beta.grid), 2)) )
  }
}
FindBeta <- data.frame(1:ngrid, RankCorr=TrueRankCorr, MuDiscr=MuDiscr)

FBeta1 <- subset(FindBeta, RankCorr > 0.80 & RankCorr < 0.83)
FBeta1$obs.num <- 1:nrow(FBeta1)
Beta.grid[FBeta1[183,1],]
Beta.grid[FBeta1[53,1],]

FBeta2 <- subset(FindBeta, RankCorr > 0.67 & RankCorr < 0.73)
FBeta2$obs.num <- 1:nrow(FBeta2)
FBeta2 <- FBeta2[order(-FBeta2$RankCorr),]
Beta.grid[FBeta2[9801,1],]
#Beta.grid[FBeta2[498,1],]

FBeta3 <- subset(FindBeta, RankCorr > 0.47 & RankCorr < 0.53)
FBeta3$obs.num <- 1:nrow(FBeta3)
FBeta3 <- FBeta3[order(-FBeta3$RankCorr),]
Beta.grid[FBeta3[45562,1],]
#Beta.grid[FBeta3[420,1],]

FBeta4 <- subset(FindBeta, RankCorr > 0.25 & RankCorr < 0.29)
FBeta4 <- FBeta4[order(-FBeta4$RankCorr),]
FBeta4$obs.num <- 1:nrow(FBeta4)
Beta.grid[FBeta4[39275,1],]
#Beta.grid[FBeta4[366,1],]

FBeta5 <- subset(FindBeta, RankCorr > 0.04 & RankCorr < 0.07)
FBeta5 <- FBeta5[order(-FBeta5$RankCorr),]
FBeta5$obs.num <- 1:nrow(FBeta5)
Beta.grid[FBeta5[980,1],]
#Beta.grid[FBeta5[248,1],]


## Define the theta.vec for each simulation setting here:
## Can probably have a separate theta.vec for the interaction terms.
BetaExt <- rbind(c(2.15, -1.0, 0.9, 0.10),
                 c(0.9, -0.1, 0.5, 0.3),
                 c(3.6, -1.5, 2.5, 2.0),
                 c(0.6, 0.05, 0.05, 0.3),
                 c(4.5, -1.7, -0.75, -2.0),
                 c(0.7, 0.1, -0.3, 0.2),
                 c(3.5, -1.25, -2.0, -2.0),
                 c(0.3, 0, -0.3, 0.35),
                 c(1.5, 1.0, -0.5, 2.5),
                 c(-0.5, -1, 0.3, 0),
                 3*c(-0.5, -1, 0.3, 0))

## Define matrix which determines correlation structure between
## conventional and novel covariates:
Btrans <- cbind(c(0.4, 0, 0, 0),
                c(0.25, 0, 0.5, 0.2))

ngrid <- 1 ## number of grid points for cross-validation
nsettings <- nrow(BetaExt)
nmethods <- 7
TrueRankCorr <- MuDiscr <- matrix(NA, nreps, nsettings)
MSELM <- MSERPEN <- MSERidge <- MSERPENAIC <- matrix(NA, nreps, nsettings)
MSEStack <- MSEDTL <- MSEAngle <- matrix(NA, nreps, nsettings)
Results <- matrix(NA, nsettings, nmethods + 3)
CorMat <- matrix(NA, nsettings, 2)
TmpCorMat <- matrix(NA, nreps, 2)

set.seed(13572)
ZBCov <- matrix(0, nrow=2, ncol=4)
betaset1 <- rep(0, p)
for(j in 1:nsettings) {
  beta.ext <- BetaExt[j,]
  beta.tl <- c(beta.ext, rep(0.0, q))
  #beta.int <- BetaInt[j]
  for(h in 1:nreps) {
    Z <- matrix(rnorm(nint*p), nrow=nint, ncol=p)
    Znew <- matrix(rnorm(nint*p), nrow=nint, ncol=p)
    Z <- scale(Z)
    Znew <- scale(Znew)
    EpsB <- matrix(rnorm(nint*q), nrow=nint, ncol=q)
    EpsBNew <- matrix(rnorm(nint*q), nrow=nint, ncol=q)

    B <- Z%*%Btrans + EpsB
    Bnew <- Znew%*%Btrans + EpsBNew
    B <- scale(B)
    Bnew <- scale(Bnew)
    X <- cbind(Z, B)
    Xnew <- cbind(Znew, Bnew)

    ZBCov <- ZBCov + crossprod(B, Z)/nint

    mu.true <- as.numeric(X%*%beta.true)
    mu.true.new <- as.numeric(Xnew%*%beta.true)
    int.mean.fn <- IntRiskScoreFn(X, gamma.true)
    int.mean.fn.new <- IntRiskScoFn(Xnew, gamma.true)
    Y <- int.mean.fn + rnorm(nint, sd=5)

    mu.ext <- as.numeric(Z%*%beta.ext)
    mu.ext.new <- as.numeric(Znew%*%beta.ext)
    external.ranking <- rank(as.numeric(mu.ext))/length(as.numeric(mu.ext))
    external.ranking.new <- rank(as.numeric(mu.ext.new))/length(as.numeric(mu.ext.new))
    Xr <- cbind(X, external.ranking)
    Xr.new <- cbind(Xnew, external.ranking.new)

    TrueRankCorr[h,j] <- cor(external.ranking, mu.true, method="spearman")
    MuDiscr[h,j] <- mean((mu.true - mu.ext)^2)

  }
}
colnames(Results) <- c("Setting", "RankCorr", "MuDiscr", "Rank Pen",
                       "Least Squares", "Ridge", "DTL", "ATL",
                       "Stacked", "Rank Pen AIC")
Results[,1] <- 1:nsettings
Results[,2] <- round(colMeans(TrueRankCorr), 2)
Results[,3] <- round(colMeans(MuDiscr), 2)
Results[,4] <- colMeans(MSERPEN)
Results[,5] <- colMeans(MSELM)
Results[,6] <- colMeans(MSERidge)
Results[,7] <- colMeans(MSEDTL)
Results[,8] <- colMeans(MSEAngle)
Results[,9] <- colMeans(MSEStack)
Results[,10] <- colMeans(MSERPENAIC)

if(nint==25) {
  fname <- "SimulationResults/Simulation2b_n25_tmp.RData"
} else if(nint==50) {
  fname <- "SimulationResults/Simulation2b_n50_tmp.RData"
} else if(nint==100) {
  fname <- "SimulationResults/Simulation2b_n100_tmp.RData"
}
#save(Results, TrueRankCorr, MuDiscr, MSERPEN, MSELM, MSERidge, MSEDTL,
#     MSEAngle, MSEStack, MSERPENAIC, file=fname)


round(Results[,1:3], 3)
#round(Results[,c(4, 10)], 3)

as.numeric(betaset1)
