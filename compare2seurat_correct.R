source('scDataAnalysis_Utilities.R')

mtx1 = readRDS('data/filtered_matrix/PBMC10k_filtered_mat.rds') ## run part3 of prepData.R to generate it
args = commandArgs(T)

frac_veg = as.numeric(args[1]) 

cost_time1 = system.time(seurat.tmp <- doBasicSeurat_new(mtx1, top_variable_features = frac_veg,
                                                        npc = 50,
                                                        reg.var = 'nCount_ATAC'))

cost_time2 = system.time(seurat.obj1 <- doBasicSeurat(mtx1, top_variable_features = frac_veg,
                                                         npc = 50,
                                                         reg.var = 'nCount_ATAC'))


## plot
seurat.tmp = FindNeighbors(seurat.tmp, reduction = 'pca', dims = 1:30,
                           k.param = 50)
seurat.tmp = FindClusters(seurat.tmp, resoltution = 0.2)

# change direction of pc if not positive correlated
pcs.tmp = seurat.tmp@reductions$pca@cell.embeddings
pcs.obj1 = seurat.obj1@reductions$pca@cell.embeddings
directions = sapply(1:ncol(pcs.obj1), function(x) 
  ifelse(cor(pcs.tmp[, x], pcs.obj1[, x]) > 0, 1, -1))
pcs.obj1 = t(t(pcs.obj1) * directions)
seurat.obj1@reductions$pca@cell.embeddings = pcs.obj1


seurat.obj1 = FindNeighbors(seurat.obj1, reduction = 'pca', dims = 1:30,
                            k.param = 50)
seurat.obj1 = FindClusters(seurat.obj1, resoltution = 0.2)


arand = adjustedRandIndex(seurat.obj1$seurat_clusters, 
                          seurat.tmp$seurat_clusters)

write(c(frac_veg, cost_time1[3], cost_time2[3], arand), 
      file = 'output/result_time_comp_withAdjRand.txt',
      append = T)



## plot results
if(T){
  dd = fread('output/result_time_comp_withAdjRand.txt')
  dd = dd[V1 %in% c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)]
  
  postscript('Figures/adjRand_compare_regressOnpca.eps', 
             width = 5, height = 5)
  
  barplot(dd$V4, ylab = 'adjust rand index', 
          names.arg = dd$V1, ylim = c(0, 1))
  dev.off()
  
  
  postscript('Figures/time_compare_regressOnpca.eps', width = 5, height = 5)
  plot(dd$V1, dd$V3, type = 'o', ylab = 'Computation time(s)', 
       col = 2, lwd = 2, xlab = 'Fraction of features used')
  lines(dd$V1, dd$V2, type = 'o', lwd = 2, col = 3)
  legend('topleft', legend = c('Seurat', 'Seurat_correct'), col = 2:3, 
         lwd = 2, bty = 'n')
  dev.off()
}
