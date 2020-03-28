library(data.table)
library(xlsx)
library(mclust)
source('scripts/dataProcessing.R')


## compare binary vs non-binary using the same clustering algorithm:louvain/kmeans
runSimpleComparison <- function(atac.mtx, 
                                nvegs = c(5000, 10000, 15000, 20000,25000), 
                                npc = 30, ctype){
  
  len = length(nvegs)
  rand.lv = matrix(0, len, 2)
  rand.kmeans = matrix(0, len, 2)
  K = length(unique(ctype))
  for(i in 1:len){
    ## non-biary
    seurat.obj <- doBasicSeurat_new(atac.mtx, npc = npc,
                                    top_variable_features = nvegs[i],
                                    norm_by = 'NA')
    seurat.obj <- FindNeighbors(seurat.obj, dims = 1:npc)
    res0 <- queryResolution4Seurat(seurat.obj, k = K, reduction = 'pca', 
                                   npc = npc)
    seurat.obj <- FindClusters(seurat.obj, resolution = res0)
    
    
    ## binary and normalized by regression
    atac.mtx.bi = 1 * (atac.mtx > 0)
    seurat.obj.bi <- doBasicSeurat_new(atac.mtx.bi, npc = npc, 
                                       top_variable_features = nvegs[i],
                                       norm_by = 'NA')
    seurat.obj.bi <- FindNeighbors(seurat.obj.bi, dims = 1:npc)
    
    res0 <- queryResolution4Seurat(seurat.obj.bi, k = K, reduction = 'pca', npc = npc)
    seurat.obj.bi <- FindClusters(seurat.obj.bi, resolution = res0)
    
    
    
    rand.lv[i, 1] = adjustedRandIndex(ctype, seurat.obj$seurat_clusters)
    rand.lv[i, 2] = adjustedRandIndex(ctype, seurat.obj.bi$seurat_clusters)
    
    
    ## comparing using kmeans 
    kmean.label <- generalCluster(seurat.obj@reductions$pca@cell.embeddings,
                                  k = K)
    kmean.label.bi <- generalCluster(seurat.obj.bi@reductions$pca@cell.embeddings,
                                     k = K)
    
    rand.kmeans[i, 1] = adjustedRandIndex(ctype, kmean.label)
    rand.kmeans[i, 2] = adjustedRandIndex(ctype, kmean.label.bi)
    
  }
  
  colnames(rand.lv) = colnames(rand.kmeans) = c('non-bi', 'bi')
  rownames(rand.lv) = rownames(rand.kmeans) = nvegs
  return(list('rand.lv' = rand.lv, 'rand.kmeans' = rand.kmeans))
}


## part 1: Generating and using synthetic data from bulk matrix (GSE74912) with noisy ####
## comparison using different number of variable genes (Fig S1B)

mtx.noisy2 <- readRDS('output/intermediate/GSE74912_noisy_q2/GSE74912_noisy_q2.rds')
ctype = sapply(colnames(mtx.noisy2), function(x) unlist(strsplit(x, '_'))[1])
rand.noisy2  = runSimpleComparison(mtx.noisy2, npc = 30, ctype = ctype)

mtx.noisy4 <- readRDS('output/intermediate/GSE74912_noisy_q4/GSE74912_noisy_q4.rds')
ctype = sapply(colnames(mtx.noisy4), function(x) unlist(strsplit(x, '_'))[1])
rand.noisy4  = runSimpleComparison(mtx.noisy4, npc = 30, ctype = ctype)

mtx.clean <- readRDS('output/intermediate/GSE74912_clean/GSE74912_clean.rds.rds')
ctype = sapply(colnames(mtx.clean), function(x) unlist(strsplit(x, '_'))[1])
rand.clean  = runSimpleComparison(mtx.clean, npc = 30, ctype = ctype)


saveRDS(rand.clean, file = 'output/intermediate/GSE74912_clean/eval_binary_rand_GSE74912_clean_PC30.rds')
saveRDS(rand.noisy2, file = 'output/intermediate/GSE74912_noisy_q2/eval_binary_rand_GSE74912_noisy_q2_PC30.rds')
saveRDS(rand.noisy4, file = 'output/intermediate/GSE74912_noisy_q4/eval_binary_rand_GSE74912_noisy_q4_PC30.rds')

nvegs = c(5000, 10000, 15000, 20000, 25000)
postscript('output/figures/GSE74912_noisy_q4/randIndex_compare_bi_nonbi_louvain_GSE74912_noisy_q4.eps',
           height = 6, width = 6)
plot(nvegs, rand.noisy4$rand.lv[, 1], pch = 17, ylim = c(0, 1),
     ylab = 'Adjusted Rand', xlab = '#variable features',
     col = 2)
points(nvegs, rand.noisy4$rand.lv[, 2], pch = 19, col = 1)
legend(x = 5000, y = 1, legend = c('non-binary', 'binary'),
       col = c(2, 1), pch = c(17, 19), bty = 'o')
dev.off()

postscript('output/figures/GSE74912_noisy_q4/randIndex_compare_bi_nonbi_kmeans_GSE74912_noisy_q4.eps',
           height = 6, width = 6)
plot(nvegs, rand.noisy4$rand.kmeans[, 1], pch = 17, ylim = c(0, 1),
     ylab = 'Adjusted Rand', xlab = '#variable features',
     col = 2)
points(nvegs, rand.noisy4$rand.kmeans[, 2], pch = 19, col = 1)
legend(x = 5000, y = 1, legend = c('non-binary', 'binary'),
       col = c(2, 1), pch = c(17, 19), bty = 'o')
dev.off()

## run all method with default setting ####
## use runMethodDefaultSetting.R to generate Fig S1A
