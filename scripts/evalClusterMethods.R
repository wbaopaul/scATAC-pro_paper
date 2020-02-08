# for clustering and/or validation
library(cluster)
library(factoextra)
library(fpc)
library(Seurat)
source('scripts/implementClusterMethods.R')

# summary each method with their default setting
sumMethods4default <- function(true.label, dtype = 'simuData', 
                               methods = c('seurat', 'monocle3', 'sc3', 
                                           'cisTopic', 'scABC', 'chromVAR', 'LSA', 'SCRAT')){
  pred.labels = list()
  rands = list()
  k = length(unique(true.label))
  for(method0 in methods){
    obj = readRDS(paste0('data/intermediate/', dtype, '/res4', method0, '_', dtype, '.rds'))
    if(grepl(method0, pattern = 'seurat', ignore.case = T)){
      # load seurat object
      obj = clust4seurat(obj, reduction = 'pca', k = k)
      pred.labels[[method0]] = obj@active.ident
    }
    
    if(method0 == 'monocle3'){
      obj <- reduceDimension(obj, max_components = 3,
                                     reduction_method = 'UMAP',
                                     verbose = F)
      
      obj <- clusterCells(obj, method = 'louvain',
                                  louvain_iter = 1, res = 0.9*10^(-2), k = 50,
                                  verbose = F)
      
      if(length(levels(pData(obj)$Cluster)) != k){
        resl0 = queryResolution4Monocle(obj, k = k, knn = 50, min_resl = 0.5)
      
        obj <- clusterCells(obj, method = 'louvain',
                          louvain_iter = 1, res = resl0 * 10^(-2),
                          verbose = F, k = 50)
      }
      pred.labels[[method0]] = pData(obj)$Cluster
    }
    
    if(method0 == 'cisTopic'){
      sele.cisTopic <- selectModel(obj, 
                                   keepBinaryMatrix = F, keepModels = F)
      cell_topic <- t(modelMatSelection(sele.cisTopic, 'cell', 'Probability'))
      pred.labels[[method0]] <- generalCluster(cell_topic, method = 'hclust', 
                                  k = k)
      
    }
    if(method0 == 'sc3'){
      pred.labels[[method0]] = obj[, 1]
    }
    
    if(method0 == 'chromVAR'){
      tsne_results <- deviationsTsne(obj, threshold = 0.5, perplexity = 10)
      pred.labels[[method0]] = cutree(hclust(dist(tsne_results)), k = k)
      
    }
    
    if(method0 %in% c('LSI', 'scABC', 'SCRAT')){
      pred.labels[[method0]] = obj
    }
    
    rands[[method0]] = adjustedRandIndex(as.character(pred.labels[[method0]]),
                                         true.label)
  }
  return(list('rands' = rands, 'pred.labels' = pred.labels))
}

updateSummary <- function(true.label, dtype, method0, rands, pred.labels){
  k = length(unique(true.label))
  obj = readRDS(paste0('data/intermediate/', dtype, '/res4', method0, '_', dtype, '.rds'))
  if(grepl(method0, pattern = 'seurat', ignore.case = T)){
    # load seurat object
    obj = clust4seurat(obj, reduction = 'pca', k = k)
    pred.labels[[method0]] = obj@active.ident
  }
  
  if(method0 == 'monocle3'){
    obj <- reduceDimension(obj, max_components = 3,
                           reduction_method = 'tSNE',
                           verbose = F)
    
    obj <- clusterCells(obj, method = 'louvain',
                        louvain_iter = 1, res = 0.9*10^(-2), k = 50,
                        verbose = F)
    
    if(length(levels(pData(obj)$Cluster)) != k){
      resl0 = queryResolution4Monocle(obj, k = k, knn = 50)
      
      obj <- clusterCells(obj, method = 'louvain',
                          louvain_iter = 1, res = resl0*10^(-2),
                          verbose = F, k = 50)
    }
    pred.labels[[method0]] = pData(obj)$Cluster
  }
  
  if(method0 == 'cisTopic'){
    sele.cisTopic <- selectModel(obj, 
                                 keepBinaryMatrix = F, keepModels = F)
    cell_topic <- t(modelMatSelection(sele.cisTopic, 'cell', 'Probability'))
    pred.labels[[method0]] <- generalCluster(cell_topic, method = 'hclust', 
                                           k = k)
    
  }
  if(method0 == 'sc3'){
    pred.labels[[method0]] = obj[, 1]
  }
  
  if(method0 == 'chromVAR'){
   
    #coVar = deviationsCovariability(obj)
    pca_coords = doDimReduction4mat(obj@assays$data$z)[[1]]
    
    pred.labels[[method0]] = cutree(hclust(dist(pca_coords)), k = k)
   
    
  }
  
  if(method0 %in% c('LSI', 'scABC', 'SCRAT')){
    pred.labels[[method0]] = obj
  }
  
  rands[[method0]] = adjustedRandIndex(as.numeric(pred.labels[[method0]]),
                                       (true.label))
  
  return(list('rands' = rands, 'pred.labels' = pred.labels))
}


getRand4SingleMethod <- function(obj, true.label, dtype = 'simuData', 
                               method0){
  k = length(unique(true.label))
  
  if(grepl(method0, pattern = 'seurat')){
    # load seurat object
    obj = clust4seurat(obj, reduction = 'pca', k = k)
    pred.labels = obj@active.ident
  }
  
  if(method0 == 'monocle3'){
    obj <- reduceDimension(obj, max_components = 2,
                           reduction_method = 'UMAP',
                           verbose = F)
    
    obj <- clusterCells(obj, method = 'louvain',
                        louvain_iter = 1, res = 0.9*10^(-2), k = 50,
                        verbose = F)
    
    if(length(levels(pData(obj)$Cluster)) != k){
      resl0 = queryResolution4Monocle(obj, k = k, knn = 50)
      
      obj <- clusterCells(obj, method = 'louvain',
                          louvain_iter = 1, res = resl0*10^(-2),
                          verbose = F, k = 50)
    }
    pred.labels = pData(obj)$Cluster
  }
  
  if(method0 == 'cisTopic'){
    sele.cisTopic <- selectModel(obj, 
                                 keepBinaryMatrix = F, keepModels = F)
    cell_topic <- t(modelMatSelection(sele.cisTopic, 'cell', 'Probability'))
    pred.labels <- generalCluster(cell_topic, method = 'hclust', 
                                             k = k)
    
  }
  if(method0 == 'sc3'){
    pred.labels = obj[, 1]
  }
  
  if(method0 == 'chromVAR'){
    
    #coVar = deviationsCovariability(obj)
    pca_coords = doDimReduction4mat(obj@assays$data$z)[[1]]
    
    pred.labels = cutree(hclust(dist(pca_coords)), k = k)
    
    
  }
  
  if(method0 %in% c('LSI', 'scABC', 'SCRAT')){
    pred.labels = obj
  }
  
  rand0 = adjustedRandIndex(as.numeric(pred.labels),
                                       true.label)
  return(rand0)
  
}


# evaluate clustering given cell*feature matrix object
eval_cluster4mat <- function(mat, cluster.label, alt.cluster.label = NULL,
                             distMethod = 'euclidean'){
  
  Dist = dist(mat, method = distMethod)
  
  label1 = as.numeric(cluster.label)
  if(!is.null(alt.cluster.label)) {
    if(class(alt.cluster.label) == 'character') alt.cluster.label = as.factor(alt.cluster.label)
    alt.cluster.label = as.numeric(alt.cluster.label)
  }
  
  res = cluster.stats(Dist, label1, alt.cluster.label)
  si = silhouette(label1, Dist)
  p <- fviz_silhouette(si, print.summary = F) + 
    scale_fill_brewer(palette = "Dark2") +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  summary.res = list('dunn_single' = res$dunn, 'dunn_ave' = res$dunn2, 
                     'rand' = res$corrected.rand, 'vi' = res$vi,
                     'silhouette' = res$clus.avg.silwidths, 'pl' = p)
  return(summary.res)
}

##calculate rand index, given true labels and npc, and a seurat object
## return data.table recoords rand index for different conditions
calRand_gTrueLabelsAndNPC <- function(seurat.obj, true.labels, npc = 20, 
                                      clust_methods = c('gcSNN', 'hclust', 'kmeans'), 
                                      reductions = c('pca', 'tsne', 'umap'), resolution = 0.1){
  
  if(is.null(seurat.obj@reductions$pca)){
    seurat.obj <- RunPCA(seurat.obj, npcs = npc, verbose = F)
  } 
  
  seurat.obj <- RunTSNE(seurat.obj, dims = 1:npc)
  seurat.obj <- RunUMAP(seurat.obj, dims = 1:npc, verbose = F)
  
  rands = matrix(0, length(clust_methods), length(reductions))
  
  for(i in 1:nrow(rands)){
    for (j in 1:ncol(rands)){
      seurat.obj <- clust4seurat(seurat.obj, npc = npc, method = clust_methods[i], 
                                 reduction = reductions[j], resolution = resolution, k = length(unique(true.labels)),
                                 clustLabelName = paste(clust_methods[i], '_', reductions[j], '_npc', npc))
      rands[i, j] = adjustedRandIndex(seurat.obj@meta.data[[paste(clust_methods[i], '_', reductions[j], '_npc', npc)]], true.labels)
    }
  }
  rownames(rands) = clust_methods
  colnames(rands) = reductions
  rands = reshape2::melt(rands, value.name = 'rand')
  rands$npc = npc
  return(rands)
} 


## compare clustering rand index using different npc, given a method and reduction
## plot tsne/umap or no plot (plotDR)
compRand_npc_gReductionAndMethod <- function(seurat.obj, clust_method = 'gcSNN', reduction = 'pca',
                                             resolution = 0.2, k = 10, npcs = c(20, 30, 50, 100), plotDR = TRUE, ...){
  
  
  if(is.null(seurat.obj@reductions$pca)){
    seurat.obj <- RunPCA(seurat.obj, npcs = max(npcs), verbose = F)
  } else if(ncol(seurat.obj@reductions$pca@cell.embeddings) < max(npcs)){
    seurat.obj <- RunPCA(seurat.obj, npcs = max(npcs), verbose = F)
  }
  
  
  i = 0
  pp_tsne = pp_umap = list()
  for(npc0 in npcs){
    i = i + 1
    seurat.obj <- clust4seurat(seurat.obj, npc = npc0, method = clust_method, 
                               reduction = reduction, resolution = resolution, k = k,
                               clustLabelName = paste(clust_method, '_', reduction, '_npc', npc0))
    if(plotDR) {
      seurat.obj = RunTSNE(seurat.obj, dims = 1:npc0)
      pp_tsne[[i]] = DimPlot(seurat.obj, reduction = 'tsne', group.by = paste(clust_method, '_', reduction, '_npc', npc0))
      seurat.obj = RunUMAP(seurat.obj, dims = 1:npc0, verbose = F)
      pp_umap[[i]] = DimPlot(seurat.obj, reduction = 'umap', group.by = paste(clust_method, '_', reduction, '_npc', npc0))
      
    }
    
  }
  #embd.mtx = seurat.obj@reductions[[reduction]]@cell.embeddings
  
  # evaluat clustering, calculate corrected rand index
  rands = rep(0, length(npcs) - 1)
  for(i in 1:length(rands)){
    #rands[i] = eval_cluster4mat(embd.mtx, 
    #seurat.obj@meta.data[[paste(clust_method, '_', reduction, '_npc', npcs[i])]],
    #seurat.obj@meta.data[[paste(clust_method, '_', reduction, '_npc', npcs[i +1])]])$rand
    rands[i] = adjustedRandIndex(seurat.obj@meta.data[[paste(clust_method, '_', reduction, '_npc', npcs[i])]],
                                 seurat.obj@meta.data[[paste(clust_method, '_', reduction, '_npc', npcs[i +1])]])
  }
  
  set.cols = brewer.pal(n = length(rands), name = 'Dark2')
  
  
  res = list('rand' = rands)
  if(plotDR){
    res$plots_tsne = pp_tsne
    res$plots_umap = pp_umap
  }
  
  return(res)
  
}


