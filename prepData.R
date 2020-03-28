source('scripts/generateSimuData.R')
source('scripts/implementClusterMethods.R')
source('scripts/dataProcessing.R')

## part1: resample data from blood cells, bulk ATAC of GSE74912 -- Fig S2 ####
## the raw matrix was resampled from bam file and 
## can be downloaded from: https://chopri.app.box.com/s/dlqybg6agug46obiu3mhevofnq4vit4t/
## under peak_by_cell_matrices/suppl_fig2/resample10k_hsc_raw_mat.rds
## suppose you saved the data in YOUR_PATH
mtx = readRDS('YOUR_PATH/resample10k_hsc_raw_mat.rds')

mtx.filter = filterMat(mtx, min_depth = 100)

rnames = rownames(mtx.filter)
mtx.filter = mtx.filter[grepl(rnames, pattern = '^chr'), ]
cnames = colnames(mtx.filter)
cnames = sapply(cnames, function(x) gsub('-', '_', x))
names(cnames) <- NULL
colnames(mtx.filter) <- cnames
dir.create('data/filtered_matrix/')
saveRDS(mtx.filter, 'data/filtered_matrix/resample10k_hsc_filtered_mat.rds')


## part2: generate synthetic data with noise (GSE74912) -- Fig S1 ####
## simulate from the bulk matrix
meta = fread('data/raw/GSE74912/summarizedMeta.txt')  ## about all health samples
## download the GSE74912_ATACseq_ALL_Counts.txt from 
## https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE74912&format=file&file=GSE74912%5FATACseq%5FAll%5FCounts%2Etxt%2Egz
tmp_file = 'data/raw/GSE74912/GSE74912_ATACseq_All_Counts.txt.gz'
download.file('https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE74912&format=file&file=GSE74912%5FATACseq%5FAll%5FCounts%2Etxt%2Egz', 
              tmp_file)
counts = fread(tmp_file)
cn = colnames(counts)
counts0 = subset(counts, select = meta$V4)
peaks = paste(counts$Chr, counts$Start, counts$End, sep = '-')
setkey(meta, V4)  
cn = colnames(counts0)
ctype = meta[cn]$V7  
uctype = unique(ctype)  
## aggregate data into group
counts0 = as.matrix(counts0)
bulk.mtx = lapply(uctype, function(x) rowSums(counts0[, ctype == x]))
bulk.mtx = do.call('cbind', bulk.mtx)
rownames(bulk.mtx) = peaks
colnames(bulk.mtx) = uctype

mtx.clean <- simulate_scatac(bulk = bulk.mtx, n_cells = 200, 
                             which_celltypes = uctype, n_frags_per_cell = 3000, 
                             rate_noise = 0, seed = 100, shuffle = FALSE)

mtx.noisy2 <- simulate_scatac(bulk = bulk.mtx, n_cells = 200, 
                              which_celltypes = uctype, n_frags_per_cell = 3000, 
                              rate_noise = 0.2, seed = 100, shuffle = FALSE)

mtx.noisy4 <- simulate_scatac(bulk = bulk.mtx, n_cells = 200, 
                              which_celltypes = uctype, n_frags_per_cell = 3000, 
                              rate_noise = 0.4, seed = 100, shuffle = FALSE)
rn = rowSums(mtx.clean)
mtx.clean = mtx.clean[rn > 10, ]
rn = rowSums(mtx.noisy2)
mtx.noisy2 = mtx.noisy2[rn > 10, ]
rn = rowSums(mtx.noisy4)
mtx.noisy4 = mtx.noisy4[rn > 10, ]
saveRDS(mtx.clean, file = 'data/filtered_matrix/GSE74912_clean/GSE74912_clean.rds')
saveRDS(mtx.noisy2, file = 'data/filtered_matrix/GSE74912_noisy_q2/GSE74912_noisy_q2.rds')
saveRDS(mtx.noisy4, file = 'data/filtered_matrix/GSE74912_noisy_q4/GSE74912_noisy_q4.rds')

## part3: prepare real data, GSE96769 -- Fig S2D ####
atac.mtx = getRealData('HSC_GSE96769', 0.01, min_depth = 500)
cn = colnames(atac.mtx)
cn[1] = gsub('#\t', '', cn[1])
colnames(atac.mtx) = cn
meta = fread('data/raw/HSC_GSE96769/ctype.txt', header = F)
atac.mtx = atac.mtx[, cn %in% meta$V1]
setkey(meta, V1)
ctype.merge = meta[colnames(atac.mtx)]$V2


saveRDS(atac.mtx, file = 'data/filtered_matrix/HSC_GSE96769_filtered_mat.rds')
saveRDS(ctype.merge, file = 'data/filtered_matrix/HSC_GSE96769_cell_type.rds')

## prepare data for Fig S3 ####
atac.mtx <- getRealData('PBMC10k', min_depth = 5000)
saveRDS(atac.mtx, file = 'data/filtered_matrix/PBMC10k_filtered_mat.rds')

