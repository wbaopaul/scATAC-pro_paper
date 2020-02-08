## get raw or simulated data (matrix) and filtering ####
source('scripts/dataProcessing.R')
data_name = 'HSC_GSE96769'

atac.mtx = readRawData(data_name, 0.01, 1000)

saveRDS(atac.mtx, file = paste0('data/intermediate/Filtered_mat/', 
                                data_name, '_filtered_mat.rds'))






## further data processing ####



## implement each method ####
# load processed data
atac.mtx <- readRDS('data/intermediate/Filtered_mat/PBMC10k_simuData.rds')
atac.mtx = atac.mtx$mtx


## evaluation metrics ####



## plot/visulize results ####

