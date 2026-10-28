
library(ggplot2)
library(patchwork)

load("DataAnalysis/AllFittedValues.RData")

rank_external <- rank(FitValMat[,"External"])
rank_ols_fitted <- rank(FitValMat[,"OLS"])
rank_ridge_fitted <- rank(FitValMat[,"Ridge"])
rank_angle_fitted <- rank(FitValMat[,"ATL"])
rank_rasps_fitted <- rank(FitValMat[,"RASPERS"])
rank_raspk_fitted <- rank(FitValMat[,"RASPERK"])

cor(rank_external, rank_ols_fitted, method="kendall")
cor(rank_external, rank_ridge_fitted, method="kendall")
cor(rank_external, rank_raspk_fitted, method="kendall")

## 0.46
## 0.45
## 0.53

nobs <- nrow(FitValMat)
meths <- c(rep("OLS", nobs), rep("Ridge Regression", nobs), rep("ATL", nobs),
           rep("RASPER (S)", nobs), rep("RASPER (K)", nobs))
nmeth <- length(unique(meths))
rank_dat <- data.frame(x=rep(rank_external, nmeth),
                       y=c(rank_ols_fitted, rank_ridge_fitted, rank_angle_fitted,
                           rank_rasps_fitted, rank_raspk_fitted),
                       `Estimation Method` = meths, check.names=FALSE)

rank_dat$`Estimation Method` <- factor(rank_dat$`Estimation Method`,
                                       levels=c("RASPER (K)", "OLS", "Ridge Regression", "ATL"))

rankdat12 <- subset(rank_dat, `Estimation Method`=="OLS" | `Estimation Method`=="RASPER (K)")
rankdat32 <- subset(rank_dat, `Estimation Method`=="Ridge Regression" | `Estimation Method`=="RASPER (K)")

#cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#D55E00")

cbPalette <- c("#D55E00", "#56B4E9")
  
  p1 <- ggplot(rankdat12) +
    aes(x, y, color = `Estimation Method`, shape=`Estimation Method`) +
    geom_abline(intercept=0, slope=1) +
    xlab("Rank of External Model") +
    ylab("Rank of Internal Model") +
    geom_point(size=0.9) +
    scale_color_manual(values = cbPalette) +
    theme_minimal() + 
    theme(axis.title = element_text(size = 7),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 7),
      legend.key.size = unit(0.4, "cm"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
      legend.box.spacing = unit(0.1, "cm")) +
    guides(color = guide_legend(override.aes = list(size = 2)))
  
  cbPalette <- c( "#D55E00", "#009E73")
  p2 <- ggplot(rankdat32) +
    aes(x, y, color = `Estimation Method`, shape=`Estimation Method`) +
    geom_abline(intercept=0, slope=1) +
    xlab("Rank of External Model") +
    ylab("Rank of Internal Model") +
    geom_point(size=0.9) +
    scale_color_manual(values = cbPalette) +
    theme_minimal() + 
    theme(axis.title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.key.size = unit(0.4, "cm"),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
          legend.box.spacing = unit(0.1, "cm")) +
    guides(color = guide_legend(override.aes = list(size = 2)))
  
  
  pdf("TablesFigures/Figure3.pdf",width=5.5, height=4.0)
  p1 + p2
  dev.off()


