## run method with default setting

source('dsAnalysis_utilities.R')

args = commandArgs(T)
dtype = args[1] ## key word for dataname
method = args[2]  ## method

kclusters = 5
if(dtype == 'simuData') mtx = readRDS('data/intermediate/Filtered_mat/PBMC10k_simuData.rds')$mtx
if(dtype == 'simuData_hard') mtx = readRDS('data/intermediate/Filtered_mat/PBMC10k_simuData_hard.rds')$mtx
if(grepl(dtype, pattern = 'resample.+hsc')) {
    mtx = readRDS(paste0('data/intermediate/Filtered_mat/', dtype, '_filtered_mat.rds'))
    kclusters = 13
}
if(grepl(dtype, pattern = 'bonemarrow')) {
  mtx = readRDS(paste0('data/raw/GSE96771/', dtype, '.rds'))
  kclusters = 6
}
mtx = 1 * mtx

stime = Sys.time()

if(method == 'monocle3') res = run_monocle3(mtx)
if(method == 'sc3') res = run_sc3(mtx, k = kclusters)
if(method == 'cisTopic') res = run_cisTopic(mtx)
if(method == 'chromVAR'){
   res = run_chromVAR(mtx, 'BSgenome.Hsapiens.UCSC.hg19')
   #res = run_chromVAR(mtx, 'BSgenome.Mmusculus.UCSC.mm9')
  
}

if(method == 'seurat') res = doBasicSeurat_new(mtx)
if(method == 'scABC') res = run_scABC(mtx, k = kclusters)
if(method == 'LSI') res = run_LSI(mtx, k = kclusters)
if(method == 'SCRAT') res = run_scrat(mtx, k = kclusters)

etime = Sys.time()
etime - stime

write(paste(method, as.numeric(etime - stime, units = 'secs')), 
      file = paste0('data/intermediate/', 'time4', dtype, '.txt'),
      append = T)

saveRDS(res, paste0('data/intermediate/', 'res4', method, '_', dtype, '.rds'))

