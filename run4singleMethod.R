source('scripts/implementClusterMethods.R')

args = commandArgs(T)
dtype = args[1] ## simuData
method = args[2]  ## method

kclusters = 5
if(dtype == 'simuData') mtx = readRDS('data/intermediate/Filtered_mat/PBMC10k_simuData.rds')$mtx
if(dtype == 'simuData_hard') mtx = readRDS('data/intermediate/Filtered_mat/PBMC10k_simuData_hard.rds')$mtx
if(dtype == 'resample10k_hsc') {
  mtx = readRDS('data/intermediate/Filtered_mat/filtered_peak_cell_mtx_200perCellType_sample10k.rds')
  kclusters = 13
}
mtx = 1 * mtx

if(method == 'monocle3') res = run_monocle3(mtx)
if(method == 'sc3') res = run_sc3(mtx, k = kclusters)
if(method == 'cisTopic') res = run_cisTopic(mtx)
if(method == 'chromVAR') res = run_chromVAR(mtx, 'BSgenome.Hsapiens.UCSC.hg19')
if(method == 'seurat') res = doBasicSeurate(mtx)
if(method == 'scABC') res = run_scABC(mtx, k = kclusters)
if(method == 'LSI') res = run_LSI(mtx, k = kclusters)
if(method == 'SCRAT') res = run_scrat(mtx, k = kclusters)


saveRDS(res, paste0('data/intermediate/', 'res4', method, '_', dtype, '.rds'))

