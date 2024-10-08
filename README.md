# Notes about cancer datasets
There are multiple cancer datasets and accessing them can be daunting. Here we provide some useful information that may help to better navigate the platforms. There are also other resources/tools dedicated to facilitate data processing, visualization and explorations of mutiple datasets such as [cBioPortal](https://www.cbioportal.org/) and [UCSC Xena](https://xenabrowser.net/datapages/). Our group is often interested in molecular subtype for breast cancer. There is a Bioconductor package in R called [AIMS](https://bioconductor.org/packages/3.18/bioc/html/AIMS.html) that could be used for that.

## The Cancer Genome Atlas [(TCGA)](https://www.cancer.gov/ccg/research/genome-sequencing/tcga)
Clinical, biospecimen, molecular characterization, and imaging data for samples from 11,000 patients (USA) spanning [33 cancer types](https://www.cancer.gov/ccg/research/genome-sequencing/tcga/studied-cancers). Including Breast Cancer (TCGA-BRCA) consisting of 1097 patients with 104 marked as deceased in vital status. Recurrence information is not available. Depending on how you download the data, the primary-keys (or identifiers) linking the different data types might be different and some-times missing. A good description about each data type is available at [UCSC Xena TCGA](https://xenabrowser.net/datapages/?cohort=TCGA%20Breast%20Cancer%20(BRCA)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
+ [GDC data portal](https://portal.gdc.cancer.gov/): By filtering and selecting different options it is possible to download many different data types.
    + **Clinical:** Downloading the clinical data directly from the portal will produce a zip file with 5 different files (clinical.tsv, exposure.tsv, family_history.tsv, follow_up.tsv, pathology_detail.tsv). The file clinical.tsv contains **survival information**, and other fields, together with basic information such as gender, age, race and ethnicity. There are two primary-key IDs in this file: **case_id** (016caf42-4e19-4444-ab5d-6cf1e76c4afa) and	**case_submitter_id** (TCGA-AO-A128). One of our former members (Jingwei) used TCGAbioLinks to retrieve other useful clinical information published in different papers (such as **molecular subtype**). The file with extended cinical information is here: [TCGA_BRCA_clinical_jingwei_2023.tsv](./TCGA_BRCA_clinical_jingwei_2023.tsv) Make sure to clean the file, as it seems to have some duplicates.
    + **Transcriptome profile:** Some data types related to gene expression include RNA-seq and miRNA-seq but mRNA-seq didn't seem to ve available at GDC. mRNA-seq is available through [cBioPortal](https://www.cbioportal.org/)
+ [TCGA Computational Tools](https://www.cancer.gov/ccg/research/genome-sequencing/tcga): Developed by TCGA network researchers and collaborators including [cBioPortal](https://www.cbioportal.org/) and [Firehose](https://gdac.broadinstitute.org/). Some members of our group found easier to download the data using cBioPortal. In the case of breast cancer you can download the full data clicking the download button from the top left corner here: [BRCA](https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018)
+ [TCGAbioLinks](https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html): This is a Bioconductor package (R language) that facilitates data retrieval. It also allows to document and automate the process. This process will keep **case_submitter_id** (TCGA-AO-A128), but **case_id** will not be included
```
# Download data for BRCA and BLCA projects
project_ids <- stringr::str_subset(TCGAbiolinks::getGDCprojects()$project_id, 'TCGA')
project_ids <- c("TCGA-BRCA”, "TCGA-BLCA")

data <- list()

for (project_id in project_ids) {
    print(project_id)
    data[[project_id]] <- TCGAbiolinks::GDCquery_clinic(project=project_id, type='clinical')
}

# Merge into single table
# (the "disease" column identifies each original table)
data <- do.call(dplyr::bind_rows, data)

# Write to file
output_path <- ‘your_path/clinical_data_BRCA_BLCA.tsv'
readr::write_tsv(data, output_path)
```

## The Pan-Cancer Analysis of Whole Genomes [(PCAWG)](https://dcc.icgc.org/pcawg)
International collaboration to identify common patterns of mutation in more than 2,600 cancer whole genomes from the International Cancer Genome Consortium (ICGC). PCAWG provides copy number, structural variation (translocations, etc.), gene expression and methylation data for many countries including USA (overlapping with TCGA data) . Detailed information about cancer and data types is available in this [table](https://dcc.icgc.org/projects/details). Donor data refers only to tumor data, while other clinical data might be related to normal samples (from blood), thus both with same icgc_donor_id but different sample/specimen IDs. Below is a detailed table for the breast cancer projects. 

DCC (Donors with molecular data), SSM (Simple Somatic Mutation), CNSM (Copy Number), StSM (Structural), SGV (Simple Germline Variation), PEXP (Protein Expression), JCN (Exon Junctions), -A, -S stands for array-based or sequencing-based

Code |Name |Site |Country |DCC |All |SSM |CNSM |StSM |SGV |METH-A |METH-S |EXP-A |EXP-S |PEXP |miRNA-S |JCN 
--- | --- | ---  | --- | ---  | --- | ---  | --- | ---  | --- | ---  | --- | ---  | --- | --- | --- | --- 
BRCA-EU|Breast ER+ and HER2- Cancer - EU/UK|Breast|EU/UK|569|569|569|344|544|--|--|--|--|--|--|--|--
BRCA-FR|Breast Cancer - FR|Breast|France|107|107|72|72|72|72|--|--|99|--|--|--|--
BRCA-KR|Breast Cancer - Very young women|Breast|South Korea|50|50|50|--|--|50|--|--|--|50|--|--|--
BRCA-UK|Breast Triple Negative/Lobular Cancer - UK|Breast|UK|150|151|141|112|30|--|--|--|--|--|--|--|--
BRCA-US|Breast Cancer - TCGA, US|Breast|US|1093|1093|1020|1045|--|--|1013|--|529|1041|298|1026|--
### Clinical data: 
It consists of 6 different files (see below). When downloading the clinical data, if more than one sample/patient was selected the information will be related to the whole data set (not only the selected patients). **Survival** information is vailable at donor.tsv
+ donor_exposure.tsv: icgc_donor_id	project_code, submitted_donor_id, exposure_type	exposure_intensity, tobacco_smoking_history_indicator, tobacco_smoking_intensity, alcohol_history, alcohol_history_intensity
+ donor_family.tsv: icgc_donor_id, project_code, submitted_donor_id, donor_has_relative_with_cancer_history, relationship_type, relationship_type_other, relationship_sex, relationship_age, relationship_disease_icd10, relationship_disease
+ donor_therapy.tsv: icgc_donor_id, project_code, submitted_donor_id, first_therapy_type, first_therapy_therapeutic_intent, first_therapy_start_interval, first_therapy_duration, first_therapy_response, second_therapy_type, second_therapy_therapeutic_intent, second_therapy_start_interval, second_therapy_duration, second_therapy_response, other_therapy, other_therapy_response
+ sample.tsv: icgc_sample_id, project_code, submitted_sample_id, icgc_specimen_id, submitted_specimen_id, icgc_donor_id, submitted_donor_id, analyzed_sample_interval, percentage_cellularity, level_of_cellularity, study
+ specimen.tsv: icgc_specimen_id, project_code, study_specimen_involved_in, submitted_specimen_id, icgc_donor_id, submitted_donor_id, specimen_type, specimen_type_other, specimen_interval, specimen_donor_treatment_type, specimen_donor_treatment_type_other, specimen_processing, specimen_processing_other, specimen_storage, specimen_storage_other, tumour_confirmed, specimen_biobank, specimen_biobank_id, specimen_available, tumour_histological_type, tumour_grading_system, **tumour_grade**, tumour_grade_supplemental, tumour_stage_system, **tumour_stage**, tumour_stage_supplemental, digital_image_of_stained_section, percentage_cellularity, level_of_cellularity
+ donor.tsv: icgc_donor_id	project_code, study_donor_involved_in, submitted_donor_id, donor_sex, **donor_vital_status**, **disease_status_last_followup** (such as relapse), donor_relapse_type, donor_age_at_diagnosis, donor_age_at_enrollment, donor_age_at_last_followup, donor_relapse_interval, donor_diagnosis_icd10, donor_tumour_staging_system_at_diagnosis, **donor_tumour_stage_at_diagnosis**, donor_tumour_stage_at_diagnosis_supplemental, **donor_survival_time**, **donor_interval_of_last_followup**, prior_malignancy, cancer_type_prior_malignancy, cancer_history_first_degree_relative
###  Expression data
+ Array-based: icgc_donor_id, project_code, icgc_specimen_id, icgc_sample_id, submitted_sample_id, analysis_id, gene_model, **gene_id**, **normalized_expression_value**, fold_change, **platform, experimental_protocol**, normalization_algorithm, other_analysis_algorithm, raw_data_repository, raw_data_accession, reference_sample_type
### Copy number data (CNSM)
Some of the relevant columns in the file are: mutation_type, copy_number, segment_mean, chromosome, chromosome_start, chromosome_end, gene_affected.
+ Mutation type: Loss, copy neutral, copy neutral LOH, gain  
+ Segment_mean (or median): Because it is for the whole segment, it can be above or below the expected number. For instance, for a loss it can be 0 and 1; but it could also be a 2. For a neutral it could be a 2; but also a 3. For a gain it is usually 3+
#### BRCA-FR (France)
Eventhough there is a file for **copy number**, "mutation type" has been recorded as **"undetermined"** for all patients and segment_mean is missing. Since vital_status has only 2 deceased (from a total of 72) **survival analysis will not be useful**. Only two additional patients with recurrence (disease_status_last_followup).
#### BRCA-EU (European Union and UK)
The name of the dataset is "Breast ER+ and HER2-" Which can be associated to **Luminal** (and normal-like?) breast cancer subtype. **Survival info is incomplete:** Vital_status is available but information about donor_survival_time and donor_interval_of_last_followup is missing. Recurrence information is also missing (disease_status_last_followup)
+ Within the same chromosome, there will be multiple entries with the same start-end for the segment (chromosome_start, chromosome_end) but associated to a different gene (gene_affected). It looks like this for the same donor/sample and at chromosome 1:
    
mutation_type |	segment_mean | chromosome_start | chromosome_end | gene_affected
--- | --- | --- | --- | ---
loss	| 1 | 75144668	| 104573697 | ENSG00000162654
loss	| 1 | 75144668	| 104573697 | ENSG00000162692

#### BRCA-KR (South Korea)
Despite that the dataset is from very young women (25 to 35 years old), **10 out of 50 died** (3 of them after going into complete remission). 10 of 50 **relapsed** and 3 of them have a status of alive. Follow-up of 1 to 15 years.

## Kaplan-Meier plotter tool [(km-plot)](https://kmplot.com/analysis/)
This is an online tool that allows you to make multiple analysis on a huge number of cancer datasets. You will be able to filter the data to select the cohort of your interest. This tools is an excellent resource to locate datasets for a specific clinical profile. At [Gyorffy 2021](https://pubmed.ncbi.nlm.nih.gov/34527184/) you will find this [table](https://github.com/Arsuaga-Vazquez-Lab/Cancer_data_notes/blob/main/gyorffy_2021_survival_55_BRCA_sets_table1.pdf) with details about 55 different breast cancer datasets especifying sample sizes, survival (Mostly overall survival and metastasis-free survival. Vital_status is often missing), recurrence and molecular subtype.Use the dataset name to locate the data within the Gene Expression Omnibus (GEO) or the European Bioinformatics Institute (EMBL-EBI).

## Molecular Taxonomy of Breast Cancer International Consortium [(METABRIC)](https://www.cbioportal.org/study/summary?id=brca_metabric)
METABRIC dataset contains gene expression, copy number alteration (CNA), single nucleotide polymorphism (SNP), and clinical data for 2000 tumors. The publically available version of the METABRIC dataset, which includes 1961 breast tumors and 548 matched normals, can be obtained from the cBioPortal link above. In the dataset, the tumors are classified into molecular subtypes using PAM50 plus Claudine-low genomic classification. The dataset includes 209 basal-like, 218 claudin-low, 224 Her2,700 Luminal A, 475 Luminal B, and 148 normal-like tumors, according to their classification (definition of subtypes and categorization results varies depending on classification methods; see below for further discussion regarding the classification methods).

## Human Breast Cell Atlas [(HBCA and iHBCA processed data)](https://cellxgene.cziscience.com/collections/48259aa8-f168-4bf5-b797-af8e88da6637)
Details about these datasets are available at [Reed at al., 2024](https://www.nature.com/articles/s41588-024-01688-9). Raw data will, eventually, be available here: [E-MTAB-13664](http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-13664/)
+ **HBCA:** [et al., 2024] made public data from single cell RNA sequencing to compile a human breast cell atlas (HBCA) assembled from 55 donors that had undergone reduction mammoplasties or risk reduction mastectomies.
+ **iHBCA:** Data from six of the largest scRNA-seq studies of the healthy breast. The iHBCA includes both fresh and frozen tissue prepared using a range of different protocols across multiple labs totaling 2.1 million cells from 286 individuals.
  
## Downloading TCGA-BRCA RNASeq, ABSOLUTE CNV, and Mutation with Clinical Information
The program "Short Program For PURPL and TP53.R". downloads TCGA-BRCA dataset including RNA-Seq, ABSOLUTE CNV, and mutation data. The program may be edited by applying different gene symbol (Gencode 36) to the variables geneName1 or geneName2. The program also attaches clinical data including survival information and molecular subtype which may be utilized in survival analysis. 

## TCGA-BRCA Survival Analysis in R
As a sample program, "TCGA_for_Survival_Analysis.R" in this repository downloads transcriptome and the copy number variation (CNV) data and performs survival analysis. The R program downloads TCGA-BRCA (Gencode 36) RNA-Seq and the ABOSLUTE CNV data from GDC using TCGAbiolinks. Users can edit the gene name and/or subset selection in a box for survival analysis. A gene name should be in the HGNC ID to correctly identify the gene. The program displays basic graphs and runs a survival analysis on the selected gene. The survival analysis identifies the optimal cut-off point that minimizes the p-value of the Log Rank Test and draws the Kaplan-Meier (KM) curve at the point. The program checks the overall trend of the difference in survival by drawing KM curves at the lower/median/upper quartile. 
