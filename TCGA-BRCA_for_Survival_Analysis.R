############For The First Time Only; Please install the libraries###############
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
BiocManager::install("SummarizedExperiment")
BiocManager::install("DESeq2")
install.packages("tidyverse")
install.packages("survminer")
install.packages("survival")
install.packages("vioplot")
install.packages('corrplot')
install.packages("gtsummary")
#########Please restart the RStudio after installing the libraries##############
################################################################################


# Please load the libraries before every use. 
library(TCGAbiolinks)
library(tidyverse)
library(SummarizedExperiment)
library(survminer)
library(survival)
library(DESeq2)
library(vioplot)
library(corrplot)
library(gtsummary)

# Get TCGA-BRCA RNA-Seq data.
queryTP <- GDCquery(project = "TCGA-BRCA",
                          data.category = "Transcriptome Profiling",
                          experimental.strategy = "RNA-Seq",
                          workflow.type = "STAR - Counts",
                          access = "open")

# Get TCGA-BRCA Absolute Copy Number data. 
queryACN <- GDCquery(project = "TCGA-BRCA",
                           data.category = "Copy Number Variation",
                           data.type = "Gene Level Copy Number",
                           workflow.type = "ABSOLUTE LiftOver",
                           access = "open")

# Download RNA-Seq and absolute copy number data.(It will take some time for the first run)
GDCdownload(queryTP)
GDCdownload(queryACN)

# Load the RNA-Seq data as summarized experiment. (This will take a several minutes)
brcaTP <- GDCprepare(queryTP, summarizedExperiment = TRUE)
# Load the RNA-Seq matrix from the summarized experiment.
brcaTPMatrix <- assay(brcaTP, "unstranded")
# Load the metadata of RNA-Seq dataset from the summarized experiment.
brcaTPGeneMetadata <- as.data.frame(rowData(brcaTP))

# Load the clinical dataset from the summarized experiment. 
clinicalData <- as.data.frame(colData(brcaTP))

# Store the clinical data and RNA-Seq matrix.
ddsTP <- DESeqDataSetFromMatrix(countData = brcaTPMatrix,
                                 colData = clinicalData,
                                 design = ~1)

# Apply the variance stabilizing transformation to the raw RNA-Seq data.
vstTP <- vst(ddsTP, blind = FALSE)

# Create the matrix.
brcaTPvst <- assay(vstTP)
###############The full normalized RNA-Seq Data is in brcaTPvst#################


# Load the aboslute copy number data. 
# Create the summarized experiment of the ACN data. (this will take some minutes)
bracaACN <- GDCprepare(queryACN, summarizedExperiment = TRUE)
# Load the ACN metadata. 
brcaACNMetadata <- as.data.frame(rowData(bracaACN))
# Create the ACN matrix.
brcaACNMatrix <- assay(bracaACN)
######The full ACN Data is in brcaACNMatrix (Many NAs at the top rows)##########


#####################Please type the first gene name here.######################
geneName1 <- "PURPL"
################################################################################

# Get the gene from the RNA-Seq dataset. 
TP_geneOfInterest <- brcaTPvst %>%
  as.data.frame()  %>%
  rownames_to_column(var = "gene_id") %>%
  gather(key = "case_id", value = "counts", -gene_id) %>%
  left_join(., brcaTPGeneMetadata, by = "gene_id") %>%
  filter(gene_name == geneName1)

# Extract patients' IDs and RNA-Seq counts. 
tempVars <- c("case_id", "counts")
TP_geneOfInterest <- TP_geneOfInterest[tempVars]

# Change colname from "counts" to "geneName_Exp"
geneNameExp1 <- paste(geneName1, "Exp", sep = "_")
colnames(TP_geneOfInterest)[2] <- geneNameExp1

# Create the overall survival using vital status, days to last followup, and days to death.
clinicalData$overallSurvival <- ifelse(clinicalData$paper_vital_status == "Alive",
                                       clinicalData$paper_days_to_last_followup,
                                       clinicalData$paper_days_to_death)

# Exclude samples missing overall survival.
clinicalData<- clinicalData[!is.na(clinicalData$overallSurvival),]

# Convert vital status. 
clinicalData$deceased <- ifelse(clinicalData$paper_vital_status== "Alive", 0, 1)

# Merge the clinical data and the gene expression data.
TCGA_geneOfInterest <- merge(TP_geneOfInterest, clinicalData, 
                             by.x = "case_id", by.y = "barcode")

# The program runs quickly from here
# Load the copy number data corresponding to the first gene expression. 
ACN_geneOfInterest <- brcaACNMatrix %>%
  as.data.frame()  %>%
  rownames_to_column(var = "gene_id") %>%
  gather(key = "case_id", value = "counts", -gene_id) %>%
  left_join(., brcaACNMetadata, by = "gene_id") %>%
  filter(gene_name == geneName1)

# Change the column names to specify the gene name. 
geneNameACN1 <- paste(geneName1, "ACN", sep = "_")
colnames(ACN_geneOfInterest) <- c("gene_id","case_id",geneNameACN1,"gene_name")

# Extract the case id and the copy number from the absolute copy number data.
myvars <- c("case_id", geneNameACN1)
ACN_geneOfInterest <- ACN_geneOfInterest[myvars]

# Truncate the case IDs to match the gene expression and the copy number data. 
TCGA_geneOfInterest$case_id <- gsub('-01.*', "", TCGA_geneOfInterest$case_id)
ACN_geneOfInterest$case_id <- gsub('-01.*', "", ACN_geneOfInterest$case_id)

# Merge the gene expression and the copy number dataset. 
TCGA_geneOfInterest <- merge(ACN_geneOfInterest, TCGA_geneOfInterest, 
                             by.x = "case_id", by.y = "case_id" )

#################The gene is prepared for Survival Analysis#####################

# Prepare a table to see the ACN for each molecular subtype.
ACNSubtypeData <- cbind.data.frame(TCGA_geneOfInterest[[geneNameACN1]],
                                   TCGA_geneOfInterest[["paper_BRCA_Subtype_PAM50"]])
colnames(ACNSubtypeData) <- c(geneNameACN1, "paper_BRCA_Subtype_PAM50")

# Display the table after running the next code. 
# Please wait some seconds to display.
ACNSubtypeData %>%
  tbl_summary(by = "paper_BRCA_Subtype_PAM50") %>%
  modify_header(label ~ "**PAM50**")


# Prepare a Pearson residual table.
ACNSubtypeData <- table(ACNSubtypeData)
chisq1 <- chisq.test(ACNSubtypeData, simulate.p.value = TRUE)
# Display the table.
corrplot(chisq1$residuals, is.cor = FALSE, addCoef.col = 'purple')

# The displays the simulated p-value of the Chi-squared in console.
chisq1$p.value

#########################Please choose the subtype.#############################
# Choose one or more from "LumA" "LumB" "Her2" "Basal" and/or "normal"
subtype <- c("LumA")
################################################################################

# Create the subset with only the subtype(s) chosen.
TCGA_Subtype <- TCGA_geneOfInterest[
  TCGA_geneOfInterest$paper_BRCA_Subtype_PAM50 %in% c(subtype),]
# Introduce NA to overall survival for removal (Warning will be displayed)
TCGA_Subtype$overallSurvival <- as.numeric(TCGA_Subtype$overallSurvival)
# Remove NA to avoid error.
TCGA_Subtype <- TCGA_Subtype[!is.na(TCGA_Subtype$overallSurvival),]

# Create labels for the violin plots
ACNLabel <- paste(geneName1, "Absolute Copy Number")
ExpLabel <- paste(geneName1, "Expression")
subtype <- paste(subtype, collapse = " ")
vioplotTitle <- paste("TCGA-BRCA", geneName1, subtype, "Gene Expression for Each ACN")
# Create the violin plot showing the gene expression distribution for each ACN
vioplot(TCGA_Subtype[[geneNameExp1]] ~ TCGA_Subtype[[geneNameACN1]], 
        col = c("#bef7ff", "#a6e2ff", "#8eccff",  "#75b7ff", "#5da1ff", "#458cff"),
        xlab = ACNLabel, ylab = ExpLabel, main = vioplotTitle)

# Run a Kendall's correlation test (ACN as the ordered factor)
# tau is the correlation coefficient (monotonic relationship), p-value included.
cor.test(TCGA_Subtype[[geneNameACN1]], TCGA_Subtype[[geneNameExp1]], method = "kendall")

###############Semi-automated Gene Expression KM-Plotter########################
#######Please enter the minimum and maximum percentile to scan cutoff###########
percentileMin <- 0.10
percentileMax <- 0.90
################################################################################

# The program scans the percentile range to find the cutoff that minimizes p-value.
# Please run the following codes. The KM-curve appears at the last code of this box.
scanBound <- c(percentileMin, percentileMax)
percentile <- quantile(TCGA_Subtype[[geneNameExp1]],scanBound)
minIndex <- percentile[[1]]
ExpStrata <- paste(geneName1, "strata", sep = "_")
TCGA_Subtype[[ExpStrata]] <- ifelse(TCGA_Subtype[[geneNameExp1]] 
                                    <= minIndex, "Low", "High")
logRank <- survdiff(Surv(TCGA_Subtype[["overallSurvival"]], 
                         TCGA_Subtype[["deceased"]]) ~ TCGA_Subtype[[ExpStrata]])
minValue <- logRank[[6]]
for(index in seq(from= percentile[[1]], to= percentile[[2]], by= 0.001))
{
  median_value <- index
  TCGA_Subtype[[ExpStrata]] <- ifelse(TCGA_Subtype[[geneNameExp1]] 
                                      <= median_value, "Low", "High")
  logRank <- survdiff(Surv(TCGA_Subtype[["overallSurvival"]], 
                           TCGA_Subtype[["deceased"]]) ~ TCGA_Subtype[[ExpStrata]])
  if(logRank[[6]] <= minValue)
  {
    minValue <- logRank[[6]]
    minIndex <- index
  }
}
TCGA_Subtype[[ExpStrata]] <- ifelse(TCGA_Subtype[[geneNameExp1]] <= minIndex, 
                                    "Low", "High")
KMTitle <- paste("TCGA-BRCA", geneName1, subtype, "Gene Expression KM Curves")
Expression <- TCGA_Subtype
fit <- survfit(Surv(Expression[["overallSurvival"]], Expression[["deceased"]]) 
               ~ Expression[[ExpStrata]] )
ggsurvplot(fit, data = Expression,
           xlab = "Days", title = KMTitle,
           pval = T,
           risk.table = T,
           tables.y.text = F)
################################################################################

# Identify each qartile of the gene expression
q <- quantile(TCGA_Subtype[[geneNameExp1]]) 
lowerQ <- q[[2]]
medianQ <- q[[3]]
higherQ <- q[[4]]

# Draw the KM-curve at the lower quartile gene expression cutoff.
TCGA_Subtype[[ExpStrata]] <- ifelse(TCGA_Subtype[[geneNameExp1]] <= lowerQ, 
                                    "Low", "High")
KMTitle <- paste("TCGA-BRCA", geneName1, subtype, 
                 "Gene Expression KM Curves at the Lower Quartile")
Expression <- TCGA_Subtype
fit <- survfit(Surv(Expression[["overallSurvival"]], Expression[["deceased"]]) ~ Expression[[ExpStrata]] )
ggsurvplot(fit, data = Expression,
           xlab = "Days", title = KMTitle,
           pval = T,
           risk.table = T,
           tables.y.text = F)


# Draw the KM-curve at the median gene expression cutoff.
TCGA_Subtype[[ExpStrata]] <- ifelse(TCGA_Subtype[[geneNameExp1]] <= medianQ, 
                                    "Low", "High")
KMTitle <- paste("TCGA-BRCA", geneName1, subtype, 
                 "Gene Expression KM Curves at the Median")
Expression <- TCGA_Subtype
fit <- survfit(Surv(Expression[["overallSurvival"]], Expression[["deceased"]]) ~ Expression[[ExpStrata]] )
ggsurvplot(fit, data = Expression,
           xlab = "Days", title = KMTitle,
           pval = T,
           risk.table = T,
           tables.y.text = F)

# Draw the KM plot at the higher quartile gene expression cutoff.
TCGA_Subtype[[ExpStrata]] <- ifelse(TCGA_Subtype[[geneNameExp1]] <= higherQ,
                                    "Low", "High")
KMTitle <- paste("TCGA-BRCA", geneName1, subtype, 
                 "Gene Expression KM Curves at the Lower Quartile")
Expression <- TCGA_Subtype
fit <- survfit(Surv(Expression[["overallSurvival"]], Expression[["deceased"]]) ~ Expression[[ExpStrata]] )
ggsurvplot(fit, data = Expression,
           xlab = "Days", title = KMTitle,
           pval = T,
           risk.table = T,
           tables.y.text = F)

# Group ACN into low ("0 or 1"), normal ("2"), higher ("3"), and very high ("4 and above")
ACN_Group <- paste(geneName1, "ACN_Group", sep = "_")
TCGA_Subtype[[ACN_Group]] <- ifelse(TCGA_Subtype[[geneNameACN1]] %in% c("1", "0"), "0 or 1", 
                                    ifelse(TCGA_Subtype[[geneNameACN1]] %in% "2", "2", 
                                           ifelse(TCGA_Subtype[[geneNameACN1]] %in% "3", "3", "4 and above")))
# Draw the KM curve with the absolute copy number
KMTitle <- paste("TCGA-BRCA", geneName1, subtype, "ACN KM Curves")
ACN <- TCGA_Subtype
fit <- survfit(Surv(ACN[["overallSurvival"]], ACN[["deceased"]]) ~ ACN[[ACN_Group]])
ggsurvplot(fit, data = ACN,
           xlab = "Days", title = KMTitle,
           pval = T,
           risk.table = T,
           tables.y.text = F)

