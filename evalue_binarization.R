library(data.table)
library(xlsx)
library(mclust)
source('scripts/dataProcessing.R')
source('dsAnalysis_utilities.R')

## part 1: use the same clustering algorithm on real data ####
## use HSC_GSE96769 data #### 
atac.mtx = readRawData('HSC_GSE96769', 0.01, min_depth = 500)
cn = colnames(atac.mtx)


## use the processed cells from the original study  ####
#orig.cellTable <- read.xlsx('data/raw/HSC_GSE96769/SupplementaryDataTable1.xlsx', 
#                            sheetName = 'PCA-TFscores.csv')
orig.cells = fread('data/raw/HSC_GSE96769/SupplementaryDataTable1.csv',
                   select = 1)$cellname
sele.cells = intersect(orig.cells, cn)
atac.mtx = atac.mtx[, sele.cells]

simp.names <- sele.cells %>% gsub("singles-", '', .) %>%
  gsub("frozen-", '', . , ignore.case = T) %>%
  gsub("scATAC-", '', ., ignore.case = T) %>%
  gsub("fresh-", '',., ignore.case = T) %>%
  gsub("160822-", '',., ignore.case = T) %>%
  gsub("160808-", '',., ignore.case = T) %>%
  gsub("160809-", '',., ignore.case = T) %>%
  gsub("160818-", '',., ignore.case = T) %>%
  gsub("160819-", '',., ignore.case = T) %>%
  gsub("20160726-", '',., ignore.case = T) %>%
  gsub("20160717-", '',., ignore.case = T) %>%
  gsub("20160617-", '',., ignore.case = T)

key1 = sapply(simp.names, function(x) unlist(strsplit(x, '-'))[1])
key2 = sapply(simp.names, function(x) unlist(strsplit(x, '-'))[2])

sample.names = ifelse(grepl(key1, pattern = '^BM|^PB'), key1, 'unknown')
ctype = ifelse(grepl(key1, pattern = '^BM|^PB'), key2, key1)
## merge GMP1,2,3 to GMP
ctype.merge = ifelse(grepl(ctype, pattern = 'GMP'), 'GMP', ctype)

saveRDS(atac.mtx, file = 'data/intermediate/Filtered_mat/HSC_GSE96769_filtered_mat.rds')
saveRDS(ctype.merge, file = 'data/intermediate/HSC_GSE96769/cell_type.rds')

## do dimension reduction and clustering using binarization or non-binarization ####


nvegs = c(5000, 10000, 15000, 20000, 25000)
len = length(fracs.veg)
rand.lv = matrix(0, len, 3)
rand.kmeans = matrix(0, len, 3)
npc = 50
for(i in 1:len){
  ## non-biary
  seurat.obj <- doBasicSeurat_new(atac.mtx, npc = npc,
                                  top.variable = nvegs[i],
                                  norm_by = 'NA')
  seurat.obj <- FindNeighbors(seurat.obj, dims = 1:npc)
  res0 <- queryResolution4Seurat(seurat.obj, k = 10, reduction = 'pca', npc = npc)
  seurat.obj <- FindClusters(seurat.obj, resolution = res0)
  
  
  ## binary and normalized by regression
  atac.mtx.bi = 1 * (atac.mtx > 0)
  seurat.obj.bi <- doBasicSeurat_new(atac.mtx.bi, npc = npc, 
                                     top.variable = nvegs[i],
                                     norm_by = 'NA')
  seurat.obj.bi <- FindNeighbors(seurat.obj.bi, dims = 1:npc)
  
  res0 <- queryResolution4Seurat(seurat.obj.bi, k = 10, reduction = 'pca', npc = npc)
  seurat.obj.bi <- FindClusters(seurat.obj.bi, resolution = res0)
 
  #using tf-idf normalization 
  seurat.obj.tfidf <- doBasicSeurat_new(atac.mtx.bi, npc = npc, 
                                     top.variable = nvegs[i],
                                     norm_by = 'tf-idf', reg.var = NULL)
  seurat.obj.tfidf <- FindNeighbors(seurat.obj.tfidf, dims = 1:npc)
  
  res0 <- queryResolution4Seurat(seurat.obj.tfidf, k = 10, reduction = 'pca', npc = npc)
  seurat.obj.tfidf <- FindClusters(seurat.obj.tfidf, resolution = res0)
  
  rand.lv[i, 1] = adjustedRandIndex(ctype.merge, seurat.obj$seurat_clusters)
  rand.lv[i, 2] = adjustedRandIndex(ctype.merge, seurat.obj.bi$seurat_clusters)
  rand.lv[i, 3] = adjustedRandIndex(ctype.merge, seurat.obj.tfidf$seurat_clusters)
  
  
  ## comparing using kmeans 
  kmean.label <- generalCluster(seurat.obj@reductions$pca@cell.embeddings,
                                k = 10)
  kmean.label.bi <- generalCluster(seurat.obj.bi@reductions$pca@cell.embeddings,
                                   k = 10)
  kmean.label.tfidf <- generalCluster(seurat.obj.tfidf@reductions$pca@cell.embeddings,
                                   k = 10)
  rand.kmeans[i, 1] = adjustedRandIndex(ctype.merge, kmean.label)
  rand.kmeans[i, 2] = adjustedRandIndex(ctype.merge, kmean.label.bi)
  rand.kmeans[i, 3] = adjustedRandIndex(ctype.merge, kmean.label.tfidf)
  
}

colnames(rand.lv) = colnames(rand.kmeans) = c('non-bi', 'bi', 'bi-tfidf')
rownames(rand.lv) = rownames(rand.kmeans) = nvegs

saveRDS(rand.lv, file = paste0('data/intermediate/HSC_GSE96769/eval_binary_rand_louvainPC', npc, '.rds'))
saveRDS(rand.kmeans, file = paste0('data/intermediate/HSC_GSE96769/eval_binary_rand_kmeansPC', npc, '.rds'))

# compare non-bi vs bi
postscript('output/figures/HSC_GSE96769/randIndex_compare_bi_nonbi_louvain.eps',
           height = 6, width = 6)
plot(nvegs, rand.lv[, 1], pch = 17, ylim = c(0, 0.3),
     ylab = 'Adjusted Rand', xlab = '#variable features',
     col = 2)
points(nvegs, rand.lv[, 2], pch = 19, col = 1)
legend(x = 19000, y = 0.1, legend = c('non-binary', 'binary'),
       col = c(2, 1), pch = c(17, 19), bty = 'o')
dev.off()

postscript('output/figures/HSC_GSE96769/randIndex_compare_bi_nonbi_kmeans.eps',
           height = 6, width = 6)
plot(nvegs, rand.kmeans[, 1], pch = 17, ylim = c(0, 0.15),
     ylab = 'Adjusted Rand', xlab = '#variable features',
     col = 2)
points(nvegs, rand.kmeans[, 2], pch = 19, col = 1)
legend(x = 19000, y = 0.025, legend = c('non-binary', 'binary'),
       col = c(2, 1), pch = c(17, 19), bty = 'o')
dev.off()

rand4ggplot <- data.table('rand' = c(rand.lv[, 1], rand.lv[, 2],
                                     rand.kmeans[, 1], rand.kmeans[, 2]),
                          'nveg' = rep(nvegs, 4),
                          'cluster' = rep(c('Louvain', 'kmeans'), each = 2 * nrow(rand.lv)),
                          'binarize' = rep(rep(c('non-binary', 'binary'), each = nrow(rand.lv))), 2 )
ggplot(data = rand4ggplot, aes(x = nveg, y = rand, 
                               color = cluster)) +
  geom_point(aes(shape = binarize))

## comparing using chromVAR ####



## part 2: using BM synthetic data with multiple clustering algorithms ####
## run simplely as in the above real data
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
                                    top.variable = nvegs[i],
                                    norm_by = 'NA')
    seurat.obj <- FindNeighbors(seurat.obj, dims = 1:npc)
    res0 <- queryResolution4Seurat(seurat.obj, k = K, reduction = 'pca', npc = npc)
    seurat.obj <- FindClusters(seurat.obj, resolution = res0)
    
    
    ## binary and normalized by regression
    atac.mtx.bi = 1 * (atac.mtx > 0)
    seurat.obj.bi <- doBasicSeurat_new(atac.mtx.bi, npc = npc, 
                                       top.variable = nvegs[i],
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

bm.noisy2 <- readRDS('data/raw/GSE96771/bonemarrow_noisy_p2.rds')
ctype = sapply(colnames(bm.noisy2), function(x) unlist(strsplit(x, '_'))[1])
rand.noisy2  = runSimpleComparison(bm.noisy2, npc = 30, ctype = ctype)

bm.noisy4 <- readRDS('data/raw/GSE96771/bonemarrow_noisy_p4.rds')
ctype = sapply(colnames(bm.noisy4), function(x) unlist(strsplit(x, '_'))[1])
rand.noisy4  = runSimpleComparison(bm.noisy4, npc = 30, ctype = ctype)

bm.clean <- readRDS('data/raw/GSE96771/bonemarrow_clean.rds')
ctype = sapply(colnames(bm.clean), function(x) unlist(strsplit(x, '_'))[1])
rand.clean  = runSimpleComparison(bm.clean, npc = 30, ctype = ctype)


saveRDS(rand.clean, file = 'data/intermediate/bonemarrow_clean/eval_binary_rand_bonemarrow_clean_PC30.rds')
saveRDS(rand.noisy2, file = 'data/intermediate/bonemarrow_noisy_p2/eval_binary_rand_bonemarrow_noisy_p2_PC30.rds')
saveRDS(rand.noisy4, file = 'data/intermediate/bonemarrow_noisy_p4/eval_binary_rand_bonemarrow_noisy_p4_PC30.rds')

nvegs = c(5000, 10000, 15000, 20000, 25000)
postscript('output/figures/bonemarrow_noisy_p4/randIndex_compare_bi_nonbi_louvain_bonemarrow_noisy_p4.eps',
           height = 6, width = 6)
plot(nvegs, rand.noisy4$rand.lv[, 1], pch = 17, ylim = c(0, 1),
     ylab = 'Adjusted Rand', xlab = '#variable features',
     col = 2)
points(nvegs, rand.noisy4$rand.lv[, 2], pch = 19, col = 1)
legend(x = 19000, y = 0.4, legend = c('non-binary', 'binary'),
       col = c(2, 1), pch = c(17, 19), bty = 'o')
dev.off()

postscript('output/figures/bonemarrow_clean/randIndex_compare_bi_nonbi_kmeans_bonemarrow_clean.eps',
           height = 6, width = 6)
plot(nvegs, rand.clean$rand.kmeans[, 1], pch = 17, ylim = c(0, 1),
     ylab = 'Adjusted Rand', xlab = '#variable features',
     col = 2)
points(nvegs, rand.clean$rand.kmeans[, 2], pch = 19, col = 1)
legend(x = 19000, y = 0.4, legend = c('non-binary', 'binary'),
       col = c(2, 1), pch = c(17, 19), bty = 'o')
dev.off()


## run all method with default setting 
## use runMethod4DefaultSetting.R



## part 3: using Erythroid synthetic data with multiple clustering algorithms ####
## run simplely as in the above real data

ery.noisy2 <- readRDS('data/raw/GSE96771/erythropoiesis_noisy_p2.rds')
ctype = sapply(colnames(ery.noisy2), function(x) unlist(strsplit(x, '_'))[1])
rand.noisy2  = runSimpleComparison(ery.noisy2, npc = 30, ctype = ctype)

ery.noisy4 <- readRDS('data/raw/GSE96771/erythropoiesis_noisy_p4.rds')
ctype = sapply(colnames(ery.noisy4), function(x) unlist(strsplit(x, '_'))[1])
rand.noisy4  = runSimpleComparison(ery.noisy4, npc = 30, ctype = ctype)

bm.clean <- readRDS('data/raw/GSE96771/erythropoiesis_clean.rds')
ctype = sapply(colnames(ery.clean), function(x) unlist(strsplit(x, '_'))[1])
rand.clean  = runSimpleComparison(ery.clean, npc = 30, ctype = ctype)


saveRDS(rand.clean, file = 'data/intermediate/erythropoiesis_clean/eval_binary_rand_erythropoiesis_clean_PC30.rds')
saveRDS(rand.noisy2, file = 'data/intermediate/erythropoiesis_noisy_p2/eval_binary_rand_erythropoiesis_noisy_p2_PC30.rds')
saveRDS(rand.noisy4, file = 'data/intermediate/erythropoiesis_noisy_p4/eval_binary_rand_erythropoiesis_noisy_p4_PC30.rds')

nvegs = c(5000, 10000, 15000, 20000, 25000)
postscript('output/figures/erythropoiesis_clean/randIndex_compare_bi_nonbi_louvain_erythropoiesis_clean.eps',
           height = 6, width = 6)
plot(nvegs, rand.clean$rand.lv[, 1], pch = 17, ylim = c(0, 1),
     ylab = 'Adjusted Rand', xlab = '#variable features',
     col = 2)
points(nvegs, rand.clean$rand.lv[, 2], pch = 19, col = 1)
legend(x = 19000, y = 0.4, legend = c('non-binary', 'binary'),
       col = c(2, 1), pch = c(17, 19), bty = 'o')
dev.off()


## run all method with default setting 
## use runMethodDefaultSetting.R

