setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/RankDataIntegration/rasper_reproduce")

load("DataAnalysis/EstimatedRegressionCoefficients.RData")

## Add target of shrinkage to regression coefficients.
beta.orig <- -c(log(2.09), log(1.63), log(1.91), log(1.57))
beta.tl <- c(beta.orig, rep(0, 6))

RegCoefs <- cbind(beta.tl, RegCoefs)
colnames(RegCoefs)[1] <- "External"
## Remove column that has augmented regression approach
RegCoefs <- RegCoefs[,-which(colnames(RegCoefs)=="Stack")]
rownames(RegCoefs) <- c("PSA $> 30$ ng/mL",
                        "Visceral Metastasis Yes",
                        "ECOG >= 2",
                        "Prior Chemo TTP",
                        "MSI Unknown",
                        "MSI Unstable",
                        "MMRd",
                        "TMB",
                        "TSG",
                        "\textit{CDK12}")

Tab4 <- RegCoefs

fname <- "TablesFigures/Table4.RData"
save(Tab4, file=fname)

library(xtable)
cap <- "Estimated regression coefficients from 3-year RMST risk models 
        for prostate cancer patients receiving an immune checkpoint 
        inhibitor after first-line chemotherapy in metastatic 
        castration-resistant prostate cancer (mCRPC) using 
        outcomes from the MSK-CHORD data."
xtab <- xtable(Tab4, digits=c(rep(1, 8)), caption=cap)

print(xtab)


