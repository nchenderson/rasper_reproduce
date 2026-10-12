################################################################
### Code to read in MSK-CHORD data from its github data source
### and clean it to create a tabular data set for prostate cancer
### patients treated with ICI.
### The cleaned data set will be saved as MSKProstateICI.csv
#################################################################

library(dplyr)
library(survival)

### 1. Read in PD-L1 Data
#PDL1 <- read.delim("~/Downloads/msk_chord_2024/data_timeline_pdl1.txt")
PDL1 <- read.delim("https://media.githubusercontent.com/media/cBioPortal/datahub/refs/heads/master/public/msk_chord_2024/data_timeline_pdl1.txt",
                   header = TRUE, sep="\t", row.names=NULL)

PDL1 <- PDL1 %>%
  select(PATIENT_ID, START_DATE, PDL1_POSITIVE)
names(PDL1) <- c("PATIENT_ID", "START_DATE_PDL1", "PDL1_POSITIVE")

## 2. Read in Mismatch Repair (MMR) Deficiency Data

#MMR <- read.delim("~/Downloads/msk_chord_2024/data_timeline_mmr.txt")
MMR <- read.delim("https://media.githubusercontent.com/media/cBioPortal/datahub/refs/heads/master/public/msk_chord_2024/data_timeline_mmr.txt",
                   header = TRUE, sep="\t", row.names=NULL)
MMR <- MMR %>%
  select(PATIENT_ID, START_DATE, MMR_ABSENT)
names(MMR) <- c("PATIENT_ID", "START_DATE_MMR", "MMR_ABSENT")

### 3. Read in genetic mutation information
#################################################
##. Note !!!
## The Mutation file has been updated a few times on github, so may not always work.
## Use a previously downloaded version on github if URL is broken
###################################################
#Mutation <- read.delim("~/Downloads/msk_chord_2024/data_mutations.txt")
Mutation <- read.delim("https://media.githubusercontent.com/media/cBioPortal/datahub/refs/heads/master/public/msk_chord_2024/data_mutation.txt",
                       header = TRUE, sep="\t", row.names=NULL)
## What is the problem with this url?
Mutation$PATIENT_ID <- substr(Mutation$Tumor_Sample_Barcode, 1, 9)

Mutation$TSG <- Mutation$Hugo_Symbol=="TP53" | Mutation$Hugo_Symbol=="PTEN" | Mutation$Hugo_Symbol=="RB1"
Mutation$CDK12ind <- Mutation$Hugo_Symbol=="CDK12"
Mutation <- Mutation %>%
  select(Hugo_Symbol, TSG, CDK12ind, Variant_Type, Tumor_Sample_Barcode, PATIENT_ID)

MutationSmall <- Mutation %>%
  group_by(PATIENT_ID) %>%
  summarize(numTSG=sum(TSG), numCDK12=sum(CDK12ind))
MutationSmall$CDK12 <- ifelse(MutationSmall$numCDK12 > 0, 1, 0)
MutationSmall$TSG <- ifelse(MutationSmall$numTSG > 0, 1, 0)

### 4. Read in ECOG scores
ecog <- read.delim("https://media.githubusercontent.com/media/cBioPortal/datahub/refs/heads/master/public/msk_chord_2024/data_timeline_performance_status.txt",
                   header = TRUE, sep="\t", row.names=NULL)
ecog <- ecog %>%
  select(PATIENT_ID, START_DATE, ECOG)
names(ecog) <- c("PATIENT_ID", "START_DATE_ECOG", "ECOG")

### 5. Read in PSA Values:
#PSAVals <- read.delim("~/Downloads/msk_chord_2024/data_timeline_psa_labs.txt")
PSAVals <- read.delim("https://media.githubusercontent.com/media/cBioPortal/datahub/refs/heads/master/public/msk_chord_2024/data_timeline_psa_labs.txt",
                      header = TRUE, sep="\t", row.names=NULL)
PSAVals <- PSAVals %>%
  select(PATIENT_ID, START_DATE, RESULT)
names(PSAVals) <- c("PATIENT_ID", "START_DATE_PSA", "PSA_VALUE")

### 6. Read in date of diagnosis information:
#DiagnosisFull <- read.delim2("~/Downloads/msk_chord_2024/data_timeline_diagnosis.txt")
DiagnosisFull <- read.delim("https://media.githubusercontent.com/media/cBioPortal/datahub/refs/heads/master/public/msk_chord_2024/data_timeline_diagnosis.txt",
                            header = TRUE, sep="\t", row.names=NULL)
DiagnosisFull <- DiagnosisFull %>%
  arrange(PATIENT_ID, START_DATE) %>%
  select(PATIENT_ID, START_DATE)
DiagnosisFirst <- DiagnosisFull[!duplicated(DiagnosisFull$PATIENT_ID),]
names(DiagnosisFirst) <- c("PATIENT_ID", "START_DATE_DIAG")

### 7. Read in clinical outcomes (OS data) data:
#ClinOutcomes <- read.delim("~/Downloads/msk_chord_2024/data_clinical_patient.txt",
#                           comment.char="#")
ClinOutcomes <- read.delim("https://media.githubusercontent.com/media/cBioPortal/datahub/refs/heads/master/public/msk_chord_2024/data_clinical_patient.txt",
                           comment.char="#")
## Create a variable indicating if a patient had visceral metastasis:
num_visc <- (ClinOutcomes$CNS_BRAIN == "Yes") + (ClinOutcomes$PLEURA == "Yes") + (ClinOutcomes$LIVER == "Yes") + (ClinOutcomes$LUNG == "Yes")
visc <- ClinOutcomes$CNS_BRAIN == "Yes" | ClinOutcomes$PLEURA == "Yes" | ClinOutcomes$LIVER == "Yes" | ClinOutcomes$LUNG == "Yes"
ClinOutcomes$VISCERAL <- ifelse(num_visc > 1, 1, 0)

ClinOutcomes <- ClinOutcomes %>%
  select(PATIENT_ID, OS_MONTHS, OS_STATUS, VISCERAL)
ClinOutcomes$OS_STATUS_NUM <- ifelse(ClinOutcomes$OS_STATUS=="1:DECEASED", 1, 0)

#### 8. Read in Sample information
#Samples <- read.delim("~/Downloads/msk_chord_2024/data_clinical_sample.txt", comment.char="#")
Samples <- read.delim("https://media.githubusercontent.com/media/cBioPortal/datahub/refs/heads/master/public/msk_chord_2024/data_clinical_sample.txt",
                      comment.char="#")

Samples$MSI <- rep("Unknown", nrow(Samples))
Samples$MSI[Samples$MSI_TYPE=="Instable"] <- "Unstable"
Samples$MSI[Samples$MSI_TYPE=="Stable"] <- "Stable"

Samples$MSI <- factor(Samples$MSI, levels=c("Unknown", "Unstable", "Stable"))
Samples <- Samples %>%
  select(PATIENT_ID, CANCER_TYPE, MSI_SCORE, MSI, TMB_NONSYNONYMOUS)
names(Samples) <- c("PATIENT_ID", "CANCER_TYPE", "MSI_SCORE", "MSI_STATUS", "TMB")

Samples <- Samples[!duplicated(Samples$PATIENT_ID),]

### 9. Read in Treatment data
#Treatment <- read.delim("~/Downloads/msk_chord_2024/data_timeline_treatment.txt")
Treatment <- read.delim("https://media.githubusercontent.com/media/cBioPortal/datahub/refs/heads/master/public/msk_chord_2024/data_timeline_treatment.txt",
                         header = TRUE, sep="\t", row.names=NULL)
Treatment <- Treatment %>%
  arrange(PATIENT_ID, START_DATE) %>%
  select(PATIENT_ID, START_DATE, SUBTYPE, AGENT)
names(Treatment) <- c("PATIENT_ID", "START_DATE_TRT", "SUBTYPE", "AGENT")

### 10. Now merge all files

Merge1 <- inner_join(ClinOutcomes, DiagnosisFirst, by="PATIENT_ID")
Merge2 <- inner_join(Merge1, Samples, by="PATIENT_ID")
ClinDiag <- inner_join(Merge2, Treatment, by="PATIENT_ID")

## Measure OS from time of diagnosis
ClinDiag$OS_FROM_DIAG <- ClinDiag$OS_MONTHS*30.4375 - ClinDiag$START_DATE_DIAG
## Make sure things are sorted by PATIENT_ID and then by START_DATE_TRT
ClinDiag <- ClinDiag %>%
  arrange(PATIENT_ID, START_DATE_DIAG)
## Create time of next treatment variable

## Need to create time to next treatment variable
## restricted to at least 30 days later.
fnext <- function(x) {
  ni <- length(x) - 1
  if(ni==0) {
    y <- NA
  } else {
    y <- rep(NA, ni + 1)
    for(k in 1:ni) {
      z <- x[x > x[k] + 30]
      if(length(z) > 0) {
        y[k] <- min(z)
      }
    }
  }
  return(y)
}

ClinDiag <- ClinDiag %>%
  group_by(PATIENT_ID) %>%
  mutate(NEXT_TRT_TIME = fnext(START_DATE_TRT)) %>%
  mutate(TRT_DIFF = NEXT_TRT_TIME - START_DATE_TRT) %>%
  ungroup()

##############################################################################
## Subset to Prostate Cancer patients and merge with Chemotherapy treatment timing
## information
ClinProstate <- subset(ClinDiag, CANCER_TYPE=="Prostate Cancer")
tmp_immuno <- subset(ClinProstate, SUBTYPE=="Immuno")
immuno_ids <- unique(tmp_immuno$PATIENT_ID)

ClinProstateImmuno <- ClinProstate[ClinProstate$PATIENT_ID %in% immuno_ids,]
ClinProstateImmuno$Chemo <- ifelse(ClinProstateImmuno$SUBTYPE=="Chemo", 1, 0)

FirstChemo <- subset(ClinProstateImmuno, Chemo==1)
FirstChemo <- FirstChemo[!duplicated(FirstChemo$PATIENT_ID),]
FirstChemo <- FirstChemo %>%
  select(PATIENT_ID, START_DATE_TRT, TRT_DIFF)
names(FirstChemo) <- c("PATIENT_ID", "START_DATE_CHEMO", "CHEMO_NEXT_TRT")
### 101 patients left after performing the above subsetting

### Now merge chemo information with a dataset that only contains treatment
## information at ICI start:
ImmunoOnly <- subset(ClinProstateImmuno, SUBTYPE=="Immuno")
ImmunoOnly <- ImmunoOnly[!duplicated(ImmunoOnly$PATIENT_ID),]

ClinProstateImmuno <- inner_join(ImmunoOnly, FirstChemo, by="PATIENT_ID")

###########################################################################
## We only want to look at ICI after first-line chemo
## So, keep observations where ICI was greater than FirstChemo
drop.ind <- ClinProstateImmuno$START_DATE_TRT < ClinProstateImmuno$START_DATE_CHEMO
ClinProstateImmuno <- ClinProstateImmuno[!drop.ind,]
## After subsetting this gives us a sample size of 85

## Create a variable which gives OS follow up in terms of days from start of ici
ClinProstateImmuno$OS_DAYS <- 30.4375*ClinProstateImmuno$OS_MONTHS
ClinProstateImmuno$OS_From_Trt <- ClinProstateImmuno$OS_DAYS - ClinProstateImmuno$START_DATE_TRT

## Now, only keep necessary variables
ClinProstateImmuno <- ClinProstateImmuno %>%
  select(PATIENT_ID, OS_From_Trt, OS_STATUS_NUM, START_DATE_TRT,
         AGENT, CHEMO_NEXT_TRT, MSI_STATUS, TMB, VISCERAL)

## Create variable capturing how long it took to progress on prior chemo
## It should be 1 if progression was 6 months and 0 if progression
## took one year or more
ClinProstateImmuno$chemo_dur <- 2 - pmin(ClinProstateImmuno$CHEMO_NEXT_TRT/180, 2)


## Get most recent PSA value for all patients
pt.list <- unique(ClinProstateImmuno$PATIENT_ID)
npts <- length(pt.list)

PSApriorICI <- ECOGpriorICI <- rep(NA, npts)
## Did not include PDL1 because the number of patients with PDL1 testing is too low
MMRpriorICI <- rep("No", npts)
for(k in 1:npts) {
  tmp_psa <- subset(PSAVals, PATIENT_ID==pt.list[k])
  tmp_ecog <- subset(ecog, PATIENT_ID==pt.list[k])
  tmp_mmr <- subset(MMR, PATIENT_ID==pt.list[k])
  tmp_clin <- subset(ClinProstateImmuno, PATIENT_ID==pt.list[k])

  if(nrow(tmp_mmr) > 0) {
    ind <- which.max(tmp_mmr$START_DATE_MMR)
    if(tmp_mmr$MMR_ABSENT[ind]) {
      MMRpriorICI[k] <- "Yes"
    }
    print(MMRpriorICI[k])
  }

  trt_start_time <- tmp_clin$START_DATE_TRT
  num_before_psa <- sum(tmp_psa$START_DATE_PSA <= trt_start_time)
  num_before_ecog <- sum(tmp_ecog$START_DATE_ECOG <= trt_start_time)

  if(num_before_psa > 0) {
    match_time_psa <- max(tmp_psa$START_DATE_PSA[tmp_psa$START_DATE_PSA <= trt_start_time])
    gg <- tmp_psa$PSA_VALUE[tmp_psa$START_DATE_PSA == match_time_psa]

    PSApriorICI[k] <- mean(tmp_psa$PSA_VALUE[tmp_psa$START_DATE_PSA == match_time_psa])
  }

  if(num_before_ecog > 0) {
    match_time_ecog <- max(tmp_ecog$START_DATE_ECOG[tmp_ecog$START_DATE_ECOG <= trt_start_time])
    ECOGpriorICI[k] <- mean(tmp_ecog$ECOG[tmp_ecog$START_DATE_ECOG == match_time_ecog])
  }
}
ClinProstateImmuno$psa <- rep("No", nrow(ClinProstateImmuno))
ClinProstateImmuno$psa[PSApriorICI > 30] <- "Yes"
ClinProstateImmuno$ecog <- ifelse(ECOGpriorICI >= 2, "Yes", "No")
ClinProstateImmuno$mmrd <- MMRpriorICI

## Merge in with mutation information here
ClinProstateImmunoNGS <- inner_join(ClinProstateImmuno, MutationSmall, by="PATIENT_ID")
ClinProstateImmunoNGS <- ClinProstateImmunoNGS %>%
  select(-AGENT, -numTSG, -numCDK12)

### You can verify the visceral metastasis information in the ClinOutcomes file
## by looking at the following:
#tumor_sites <- read.delim("~/Downloads/msk_chord_2024/data_timeline_tumor_sites.txt")

#pt.list <- ClinProstateImmunoNGS$PATIENT_ID
#npts <- length(pt.list)
#ViscMets <- rep(NA, npts)
#for(k in 1:npts) {
#    tst_tumor <- subset(tumor_sites, PATIENT_ID==pt.list[k])
#    nbrain <- sum(tst_tumor$TUMOR_SITE == "CNS/Brain")
#    nlung <- sum(tst_tumor$TUMOR_SITE == "Lung")
#    nliver <- sum(tst_tumor$TUMOR_SITE == "Liver")
#    npleura <- sum(tst_tumor$TUMOR_SITE == "Pleura")

#    ViscMets[k] <- ifelse(nbrain > 0 | nlung > 0 | nliver > 0 | npleura > 0, 1, 0)
#}

### Create inverse censoring weights and save results:
ici_cens <- survfit(Surv(OS_From_Trt, 1 - OS_STATUS_NUM) ~ 1, data=ClinProstateImmunoNGS)
Ghat <- approxfun(ici_cens$time, ici_cens$surv)
ICWeights <- 1/Ghat(ClinProstateImmunoNGS$OS_From_Trt)
ClinProstateImmunoNGS$ICWeights <- ICWeights

### Also save RMST pseudo-values
sfit <- survfit(Surv(OS_From_Trt, OS_STATUS_NUM) ~ 1, data=ClinProstateImmunoNGS)
## Do 1 or 2 or 3 year RMST pseudo-values.

ypseudo_3yr <- pseudo(sfit, times=365.25*3, type="RMST")
ClinProstateImmunoNGS$ypseudo_3yr <- ypseudo_3yr

## Remove patient with missing time on chemo information:
ClinProstateImmunoNGS <- ClinProstateImmunoNGS[!is.na(ClinProstateImmunoNGS$chemo_dur),]

### Now, save file as MSKProstateICI.csv in whichever local directory you prefer
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/RankDataIntegration/rasper_reproduce")
write.csv(ClinProstateImmunoNGS, file="DataAnalysis/MSKProstateICI.csv",
          row.names=FALSE)





