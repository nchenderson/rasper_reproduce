
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/RankDataIntegration/rasper_reproduce")
load("SimulationResults/Simulation1b_n50.RData")

Tab2 <- Results
OLS_MSE <- Tab2[,5]
Tab2[,4] <- Tab2[,4]/OLS_MSE
#Tab2[,5] <- colMeans(MSERPENKendall)/OLS_MSE
Tab2[,5] <- Tab2[,11]/OLS_MSE
for(k in 6:10) {
  Tab2[,k] <- Tab2[,k]/OLS_MSE
}
Tab2 <- Tab2[,-11]
colnames(Tab2)[5] <- "Rank Pen Kendall"
round(Tab2, 3)

fname <- "TablesFigures/Table2.RData"
save(Tab2, file=fname)

library(xtable)
cap <- "Simulation Study 1(b). Relative mean-squared error is reported for each method. 
        The relative MSE is defined as the MSE divided by the MSE obtained by least-squares estimation.
        For each of the 11 settings, the expected rank correlation ? and distance between the 
        internal and external mean functions muI and muE are reported."
xtab <- xtable(Tab2, digits=c(0, 0, 2, 1, rep(2, 7)), caption=cap)

print(xtab, include.rownames = FALSE)



