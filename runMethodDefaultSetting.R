## run method with default setting

source('scripts/implementClusterMethods.R')

args = commandArgs(T)
dtype = args[1] ## key word for data name: one of HSC_GSE96769, resample10k_hsc,  GSE74912_clean, GSE74912_noisy_q2, GSE74912_noisy_q2_q4
method = args[2]  ## clustering method

kclusters = 13
mtx = readRDS(paste0('data/filtered_mat/', dtype, '_filtered_mat.rds'))

if(grepl(dtype, pattern = 'GSE74912')) {
  ## systhetic data to study binary vs non-binary
  mtx = readRDS(paste0('output/intermediate/', dtype, '/', dtype, '.rds'))
}

if(grepl(dtype, pattern = 'HSC_GSE96769')) {
  kclusters = 10
}

mtx = 1 * mtx

stime = Sys.time()

if(method == 'monocle3') res = run_monocle3(mtx)
if(method == 'sc3') res = run_sc3(mtx, k = kclusters)
if(method == 'cisTopic') res = run_cisTopic(mtx)
if(method == 'chromVAR'){
   gname =  'BSgenome.Hsapiens.UCSC.hg19'
   res = run_chromVAR(mtx, gname)
}

## note to change top_variable_features when differnt # of variable genes was needed
if(method == 'seurat') res = doBasicSeurat(mtx, 5000)
if(method == 'seurat_correct') {
  res = doBasicSeurat_new(mtx, top_variable_features = 5000) 
}

if(method == 'scABC') res = run_scABC(mtx, k = kclusters)
if(method == 'LSI') res = run_LSI(mtx, k = kclusters)
if(method == 'SCRAT') res = run_scrat(mtx, k = kclusters)

etime = Sys.time()
etime - stime

## save result
dir.create(paste0('output/intermediate/', dtype), showWarnings = F)

write(paste(method, as.numeric(etime - stime, units = 'secs')), 
      file = paste0('output/intermediate/', dtype, '/time4', dtype, '.txt'),
      append = T)


saveRDS(res, paste0('output/intermediate/', dtype, '/res4', method, '_', dtype, '.rds'))

