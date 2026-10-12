setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/RankDataIntegration/rankingshrink_reproduce")
#load("SimulationResults/Simulation1b_kendall_n50.RData")
load("SimulationResults/Simulation2a_n50.RData")

Tab3 <- Results
OLS_MSE <- Tab3[,5]
Tab3[,4] <- Tab3[,4]/OLS_MSE
#Tab2[,5] <- colMeans(MSERPENKendall)/OLS_MSE
Tab3[,5] <- Tab3[,11]/OLS_MSE
for(k in 6:10) {
  Tab3[,k] <- Tab3[,k]/OLS_MSE
}
Tab3 <- Tab3[,-11]
colnames(Tab3)[5] <- "Rank Pen Kendall"
round(Tab3, 3)

fname <- "TablesFigures/Table3.RData"
save(Tab3, file=fname)

library(xtable)
cap <- "Simulation Study 2. Relative mean-squared error is reported for each method.
        The relative MSE is defined as the MSE divided by the MSE obtained by
        least-squares estimation. For each of the 10 simulation settings, the
        expected rank correlation ? and expected distance
        ? between the internal and external mean
        functions ? and ? are reported, and the relative MSE of
        the best performing method in each setting is in bold text.."
xtab <- xtable(Tab3, digits=c(0, 0, 2, 1, rep(2, 7)), caption=cap)

print(xtab, include.rownames = FALSE)
