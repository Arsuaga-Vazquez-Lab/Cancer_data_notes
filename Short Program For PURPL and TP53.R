##########For The Very First Time Only; Please install the libraries############
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("TCGAbiolinks")
# BiocManager::install("SummarizedExperiment")
# BiocManager::install("DESeq2")
# install.packages("tidyverse")
# install.packages("survminer")
# install.packages("survival")
# install.packages("vioplot")
# install.packages('corrplot')
install.packages('Hmisc')
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


# Download RNA-Seq and absolute copy number data.
GDCdownload(queryTP)
GDCdownload(queryACN)



# Load the RNA-Seq data as summarized experiment.
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

# Load the aboslute copy number data. 
# Create the summarized experiment of the copy number data and create its dataframe. 
bracaACN <- GDCprepare(queryACN, summarizedExperiment = TRUE)
brcaACNMatrix <- assay(bracaACN)
# Load the absolute copy number metadata. 
brcaACNMetadata <- as.data.frame(rowData(bracaACN))

# Get TCGA-BRCA Mutation data.
queryMut <- GDCquery(project = "TCGA-BRCA",
                     data.category = "Simple Nucleotide Variation",
                     data.type = "Masked Somatic Mutation",
                     workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
                     access = "open")

# Download the data
GDCdownload(queryMut)


brcaMut <- GDCprepare(queryMut, summarizedExperiment = TRUE)

queryMethyl <- GDCquery(project = "TCGA-BRCA",
                        data.category = "DNA Methylation",
                        data.type = "Methylation Beta Value",
                        platform = "Illumina Human Methylation 450",
                        access = "open")

# Download the data
GDCdownload(queryMethyl)

# # Prepare the data
methylation_data <- GDCprepare(queryMethyl)

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

# Create overall survival information using vital status, days to last followup, and days to death.
clinicalData$overallSurvival <- ifelse(clinicalData$paper_vital_status == "Alive",
                                       clinicalData$paper_days_to_last_followup,
                                       clinicalData$paper_days_to_death)

# Exclude samples missing overall survival.
clinicalData<- clinicalData[!is.na(clinicalData$overallSurvival),]

# Convert vital status. 
clinicalData$deceased <- ifelse(clinicalData$paper_vital_status== "Alive", 0, 1)

# Merge the clinical data and the gene expression data.
clinicalData <- clinicalData[!duplicated(clinicalData$barcode), ]
TP_geneOfInterest <- TP_geneOfInterest[!duplicated(TP_geneOfInterest$case_id), ]
TCGA_geneOfInterest <- merge(TP_geneOfInterest, clinicalData, by.x = "case_id", by.y = "barcode")

###################Please type the second gene name here.#######################
geneName2 <- "TP53"
################################################################################

# Get the gene from the RNA-Seq dataset. 
TP_geneOfInterest <- brcaTPvst %>%
  as.data.frame()  %>%
  rownames_to_column(var = "gene_id") %>%
  gather(key = "case_id", value = "counts", -gene_id) %>%
  left_join(., brcaTPGeneMetadata, by = "gene_id") %>%
  filter(gene_name == geneName2)

# Extract patients' IDs and RNA-Seq counts. 
tempVars <- c("case_id", "counts")
TP_geneOfInterest <- TP_geneOfInterest[tempVars]

# Change colname from "counts" to "geneName_Exp"
geneNameExp2 <- paste(geneName2, "Exp", sep = "_")
colnames(TP_geneOfInterest)[2] <- geneNameExp2
# Merge the clinical data and the gene expression data.
TCGA_geneOfInterest <- TCGA_geneOfInterest[!duplicated(TCGA_geneOfInterest$case_id), ]
TP_geneOfInterest <- TP_geneOfInterest[!duplicated(TP_geneOfInterest$case_id), ]
TCGA_geneOfInterest <- merge(TP_geneOfInterest, TCGA_geneOfInterest, by = "case_id")
################################################################################
# Mutation Data
################################################################################
Mut_geneOfInterest <- brcaMut %>%
  filter(Hugo_Symbol == geneName2)
MutVector <- Mut_geneOfInterest$Tumor_Sample_Barcode
# MutVector <- gsub('-01.*', "", MutVector)
# MutVector <- gsub('-06.*', "", MutVector)
# TCGA_geneOfInterest$case_id <- gsub('-01.*', "", TCGA_geneOfInterest$case_id)
# Truncate the case ids to match the case ids of the gene expression and the copy number datasets.

truncate_string <- function(input_string, delimiter) {
  match_position <- regexpr(delimiter, input_string)
  if (match_position == -1) {
    return(input_string)  # delimiter not found, return the original string
  }
  end_position <- match_position + attr(match_position, "match.length") - 1
  result <- substring(input_string, 1, end_position)
  return(result)
}

truncate_vector <- function(input_vector, delimiter) {
  truncated_vector <- sapply(input_vector, truncate_string, delimiter = delimiter)
  return(truncated_vector)
}

input_vector <- MutVector
delimiter <- "-01A"
mut_vector <- truncate_vector(input_vector, delimiter)
result_vector <- as.data.frame(mut_vector)

TCGA_geneOfInterest$TP53_Mut_Status <- ifelse(TCGA_geneOfInterest$sample %in% mut_vector, "TP53_Mut", "TP53_WT")
################################################################################
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

# Extract the case id and the copy number from the absolute copy number dataset.
myvars <- c("case_id", geneNameACN1)
ACN_geneOfInterest <- ACN_geneOfInterest[myvars]

# Truncate the case ids to match the case ids of the gene expression and the copy number datasets. 
TCGA_geneOfInterest$case_id <- gsub('-01.*', "", TCGA_geneOfInterest$case_id)
ACN_geneOfInterest$case_id <- gsub('-01.*', "", ACN_geneOfInterest$case_id)
ACN_geneOfInterest$case_id <- gsub('-06.*', "", ACN_geneOfInterest$case_id)
testVector <- as.data.frame(TCGA_geneOfInterest$case_id)
testVector <- as.data.frame(ACN_geneOfInterest$case_id)

# Merge the gene expression and the copy number dataset. 
TCGA_geneOfInterest <- TCGA_geneOfInterest[!duplicated(TCGA_geneOfInterest$case_id), ]
ACN_geneOfInterest <- ACN_geneOfInterest[!duplicated(ACN_geneOfInterest$case_id), ]
TCGA_geneOfInterest <- merge(ACN_geneOfInterest, TCGA_geneOfInterest, by = "case_id")

################################################################################

# Load the copy number data corresponding to the second gene expression. 
ACN_geneOfInterest <- brcaACNMatrix %>%
  as.data.frame()  %>%
  rownames_to_column(var = "gene_id") %>%
  gather(key = "case_id", value = "counts", -gene_id) %>%
  left_join(., brcaACNMetadata, by = "gene_id") %>%
  filter(gene_name == geneName2)

# Change the column names to specify the gene name. 
geneNameACN2 <- paste(geneName2, "ACN", sep = "_")
colnames(ACN_geneOfInterest) <- c("gene_id","case_id",geneNameACN2,"gene_name")

# Extract the case id and the copy number from the absolute copy number dataset.
myvars <- c("case_id", geneNameACN2)
ACN_geneOfInterest <- ACN_geneOfInterest[myvars]

# Truncate the case ids to match the case ids of the gene expression and the copy number datasets. 
ACN_geneOfInterest$case_id <- gsub('-01.*', "", ACN_geneOfInterest$case_id)

# Merge the gene expression and the copy number dataset. 
TCGA_geneOfInterest <- TCGA_geneOfInterest[!duplicated(TCGA_geneOfInterest$case_id), ]
ACN_geneOfInterest <- ACN_geneOfInterest[!duplicated(ACN_geneOfInterest$case_id), ]
TCGA_geneOfInterest <- merge(ACN_geneOfInterest, TCGA_geneOfInterest, by = "case_id" )
#########################Please choose the subtype.#############################
# Choose one or more from "LumA", "LumB", "Her2", "Basal", and/or "Normal"
subtype <- c("LumA")
################################################################################
TCGA_Subtype <- TCGA_geneOfInterest[TCGA_geneOfInterest$paper_BRCA_Subtype_PAM50 
                                    %in% c(subtype),]

subtype <- paste(subtype, collapse = " ")
TCGA_Subtype$overallSurvival <- as.numeric(TCGA_Subtype$overallSurvival)
TCGA_Subtype <- TCGA_Subtype[!is.na(TCGA_Subtype$overallSurvival),]