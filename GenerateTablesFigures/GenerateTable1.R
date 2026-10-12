
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/RankDataIntegration/rankingshrink_reproduce")

load("SimulationResults/Simulation1a_n50.RData")

Tab1 <- Results
OLS_MSE <- Tab1[,5]
Tab1[,4] <- Tab1[,4]/OLS_MSE
#Tab1[,5] <- colMeans(MSERPENKendall)/OLS_MSE
Tab1[,5] <- Tab1[,11]/OLS_MSE
for(k in 6:10) {
    Tab1[,k] <- Tab1[,k]/OLS_MSE
}
Tab1 <- Tab1[,-11]

colnames(Tab1)[5] <- "Rank Pen Kendall"
round(Tab1, 3)

fname <- "TablesFigures/Table1.RData"
save(Tab1, file=fname)

library(xtable)
cap <- "Simulation Study 1(a). Relative mean-squared error is reported for each method. 
        The relative MSE is defined as the MSE divided by the MSE obtained by least-squares estimation.
        For each of the 11 settings, the expected rank correlation ? and distance between the 
        internal and external mean functions muI and muE are reported."
xtab <- xtable(Tab1, digits=c(0, 0, 2, 1, rep(2, 7)), caption=cap)
print(xtab, include.rownames = FALSE)
