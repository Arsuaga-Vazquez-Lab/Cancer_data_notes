# Notes about cancer datasets
There are multiple cancer datasets and accessing them can be daunting. Here we provide some useful information that may help to better navigate the platforms.  

## The Cancer Genome Atlas [(TCGA)](https://www.cancer.gov/ccg/research/genome-sequencing/tcga)
Clinical, biospecimen, molecular characterization, and imaging data for samples from 11,000 patients (USA) spanning [33 cancer types](https://www.cancer.gov/ccg/research/genome-sequencing/tcga/studied-cancers). Including Breast Cancer (TCGA-BRCA) with 1000+ patients. Depending on how you download the data, the primary-keys (or identifiers) linking the different data types might be different and some-times missing.
+ [GDC data portal](https://portal.gdc.cancer.gov/): By filtering and selecting different options it is possible to download many different data types. Downloading the clinical data directly from the portal will produce a zip file with 5 different files (clinical.tsv, exposure.tsv, family_history.tsv, follow_up.tsv, pathology_detail.tsv). The file clinical.tsv contains **survival information**, and other fields, together with basic information such as gender, age, race and ethnicity. There are two primary-key IDs in this file: **case_id** (016caf42-4e19-4444-ab5d-6cf1e76c4afa) and	**case_submitter_id** (TCGA-AO-A128)
+ [TCGAbioLinks](https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html)
This is a Bioconductor package (R language) that facilitates data retrieval. It also allows to document and automate the process. This process will keep **case_submitter_id** (TCGA-AO-A128), but **case_id** will not be included
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
International collaboration to identify common patterns of mutation in more than 2,600 cancer whole genomes from the International Cancer Genome Consortium (ICGC). PCAWG provides copy number, structural variation (translocations, etc.), gene expression and methylation data for many countries including USA (overlapping with TCGA data) . Detailed information about cancer and data types is available in this [table](https://dcc.icgc.org/projects/details).
