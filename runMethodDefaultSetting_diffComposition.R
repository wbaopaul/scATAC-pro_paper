## run method with default setting and benchmarking on different compositions
library(gtools)

source('scripts/implementClusterMethods.R')
source('scripts/evalClusterMethods.R')
source('scripts/dataProcessing.R')
args = commandArgs(T)
dtype = args[1] ## simuData
method = args[2]  ## method


if(dtype %in% c('simuData', 'simuData_hard')) {
  mtx = readRDS(paste0('data/intermediate/Filtered_mat/PBMC10k_', dtype, '.rds'))
  true.label = mtx$cluster_labels$Group
  mtx = mtx$mtx
  cnames = colnames(mtx)
  kclusters = 5
}

if(grepl(dtype, pattern = 'resample')) {
  mtx = readRDS(paste0('data/intermediate/Filtered_mat/', dtype, '_filtered_mat.rds'))
  kclusters = 13
  cnames = colnames(mtx)
  true.label = sapply(cnames, function(x) unlist(strsplit(x, '-'))[1])
}
mtx = 1 * mtx

stime = Sys.time()
## generate different composition
set.seed(2019)
props = rdirichlet(100, rep(3, 13))
N = 1000 # total number of cells
groups = unique(true.label)
rands = rep(0, nrow(props))
for(i in 1:nrow(props)){
  # generate data
  nn = floor(N * props[i, ])
  nn[nn > 200] = 200
  nn[nn < 10] = 10
  cnames.sele = lapply(1:kclusters, function(x) {
                          cnames0 = cnames[true.label == groups[x]]; return(cnames0[1:nn[x]])})
  cnames.sele = do.call('c', cnames.sele)
  mtx.sele = mtx[, cnames %in% cnames.sele]
  mtx.sele = filterMat(mtx.sele, min_depth = 1000)
  true.label.sele = sapply(colnames(mtx.sele), function(x) unlist(strsplit(x, '-'))[1])
  
  if(method == 'monocle3') res = run_monocle3(mtx.sele)
  if(method == 'sc3') res = run_sc3(mtx.sele, k = kclusters)
  if(method == 'cisTopic') res = run_cisTopic(mtx.sele, nCores = 8)
  if(method == 'chromVAR') {
    res = run_chromVAR(mtx.sele, 'BSgenome.Hsapiens.UCSC.hg19')
  }
  if(method == 'seurat') res = doBasicSeurate(mtx.sele)
  if(method == 'seurat_correct') res = doBasicSeurat_new(mtx.sele)
  if(method == 'scABC') res = run_scABC(mtx.sele, k = kclusters)
  if(method == 'LSI') res = run_LSI(mtx.sele, k = kclusters)
  if(method == 'SCRAT') res = run_scrat(mtx.sele, k = kclusters)
  rands[i] = getRand4SingleMethod(res, method0 = method, true.label = true.label.sele)
  
}

etime = Sys.time()
etime - stime
write(paste(method, as.numeric(etime - stime, units = 'secs')), 
      file = paste0('data/intermediate/', dtype, '/time4', dtype, '_diffComposition.txt'),
      append = T)


write.table(rands, file = paste0('data/intermediate/randIndexs/', 'rand4', method, '_', dtype, '.txt'),
            row.names = F, quote = F, col.names = F, sep = '\t')

