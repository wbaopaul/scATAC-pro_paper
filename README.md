scATAC-pro_study
-----------------

Data and scripts for manuscript "scATAC-pro: a comprehensive workbench for single-cell chromatin accessibility sequencing data"


## Scripts


- prepData.R --- Script to download/generat and process simulated and real data

- scripts/   --- Includes all scripts to process data, implement and evaluate different clustering methods

- runMethodDefaultSetting.R --- Perform clustering analysis given a dataset and clustering method (Fig S1A, S2A, S2B, S2D)

- evalue_binarization.R  --- Study the effect of data binarization on clustering analysis (Fig S1B)

- runMethodDefaultSetting_differentComp.R --- Perform clustering analysis for simulated data vary in cell type composition (Fig S2C)

- summarizeRes.R --- Summarize all intermediate results and plot figure (Fig S1, S2)

- compare2seurat_correct.R  -- Compare PCA implemented in Seurat and scATAC-pro (Fig S3)


## Data

### Data for benchmarking
Run prepData.R to prepare all data for the analysis. The filtered peak-by-cell matrices are saved under data/filtered_matrix, except the PBMC data was saved [here](https://chopri.app.box.com/s/dlqybg6agug46obiu3mhevofnq4vit4t/). 

### Data for case study
The processed data for three case studies along with some annotation files for scATAC-pro were saved [here](https://chopri.app.box.com/s/dlqybg6agug46obiu3mhevofnq4vit4t/).   


## Citation

Yu W, Uzun Y, Zhu Q, Chen C, Tan K. *scATAC-pro: a comprehensive workbench for single-cell chromatin accessibility sequencing data.* bioRxiv.org; 2019
doi: https://doi.org/10.1101/824326


