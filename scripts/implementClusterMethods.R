library(data.table)
library(compiler)
library(magrittr)
library(mclust)
library(dplyr)

library(chromVAR)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Mmusculus.UCSC.mm9)
library(SummarizedExperiment)
library(BiocParallel)
library(JASPAR2016)

library(Rtsne)
library(Seurat)
library(DropletUtils)
library(cisTopic)
#library(scABC)

library(monocle)
suppressMessages(library(flexclust))
suppressMessages(library(mcclust))
#library(reticulate)
#import("louvain")

library(SC3)

library(reticulate)
use_python('/mnt/isilon/tan_lab/yuw1/local_tools/anaconda3/bin/python')

doDimReduction4mat <- function(mtx, max_pc = 20, doTSNE = F){
  ## DO SVD on tf_idf normalized matrix                             
  set.seed(1234)
  SVDtsne = irlba(mtx, max_pc, max_pc, maxit=1000)
  d_diagtsne = matrix(0, nrow=length(SVDtsne$d), ncol=length(SVDtsne$d))
  diag(d_diagtsne) = SVDtsne$d
  SVDtsne_vd = t(d_diagtsne %*% t(SVDtsne$v))
  rownames(SVDtsne_vd) = colnames(mtx)
  colnames(SVDtsne_vd) = paste0('pca_', 1:ncol(SVDtsne_vd))
  
  ## Run TSNE to 2 dimensions
  tsne_coords = NULL
  if(doTSNE){
    tsnetfidf = Rtsne(SVDtsne_vd, pca=F, perplexity=30, max_iter=5000)
    
    tsne_coords = as.data.frame(tsnetfidf$Y)
    colnames(tsne_coords) = c('tsne_1', 'tsne_2')
    rownames(tsne_coords) = colnames(mtx) 
  }
  
  pca_coords = SVDtsne_vd
  return(list('pca_coords' = pca_coords, 'tsne_coords' = tsne_coords))
}

## implement of scRNA-seq specific methods ####
# do normalization, pca using Seurat, with default setting
doBasicSeurate <- function(mtx, npc = 50, top.variable = NULL,
                           doScale = T, doCenter = T, doLog = T, assay = 'ATAC', 
                           regr.var = 'nCount_ATAC'){
  
  # top.variabl -- use top most variable features
  if(doLog) mtx = log2(mtx + 1)
  seurat.obj = CreateSeuratObject(mtx, project = 'scATAC', assay = assay,
                                  names.delim = '-', min.cells = 10, min.features = 10)
  cell.names = colnames(mtx)
  
  
  seurat.obj <- FindVariableFeatures(object = seurat.obj, 
                                     nfeatures = ifelse(is.null(top.variable), 5000, 
                                                        floor(nrow(mtx) * top.variable)))
  seurat.obj <- ScaleData(object = seurat.obj, 
                          features = VariableFeatures(seurat.obj), 
                          vars.to.regress = regr.var, do.scale = doScale, do.center = doCenter)
  
  
  seurat.obj <- RunPCA(object = seurat.obj, 
                       features = VariableFeatures(object = seurat.obj),
                       verbose = FALSE, seed.use = 10, npc = npc)
  
  
  return(seurat.obj)
}
doBasicSeurate = cmpfun(doBasicSeurate)


regress_on_pca <- function(seurat.obj, reg.var = 'nCount_ATAC'){
  
  pcs = seurat.obj@reductions$pca@cell.embeddings
  pcs.reg = pcs
  for(i in 1:length(reg.var)){
    
    reg.var0 = seurat.obj[[reg.var[i]]][[1]]
    pcs.reg = apply(pcs.reg, 2, function(x) lm(x ~ reg.var0)$residual )
    
  }
  colnames(pcs.reg) = colnames(pcs)
  seurat.obj@reductions$pca@cell.embeddings = pcs.reg
  return(seurat.obj)
}


#could be normalized by log, tf-idf or NA
doBasicSeurat_new <- function(mtx, npc = 50, top_variable_features = 0.2, 
                              doScale = T, doCenter = T, assay = 'ATAC',
                              reg.var = NULL, norm_by = 'log', project = 'scATAC'){
  
  # top.variabl -- use top most variable features
  seurat.obj = CreateSeuratObject(mtx, project = project, assay = assay,
                                  names.delim = '-')
  
  if(norm_by == 'log') seurat.obj@assays$ATAC@data <- log1p(seurat.obj@assays$ATAC@data)/log(2)
  if(norm_by == 'tf-idf') seurat.obj@assays$ATAC@data <- TF.IDF(seurat.obj@assays$ATAC@data)
  
  nveg = ifelse(top_variable_features > 1, top_variable_features, floor(nrow(mtx) * top_variable_features))
  seurat.obj <- FindVariableFeatures(object = seurat.obj,
                                     selection.method = 'vst',
                                     nfeatures = nveg)
  seurat.obj <- ScaleData(object = seurat.obj,
                          features = VariableFeatures(seurat.obj),
                          vars.to.regress = NULL, do.scale = doScale,
                          do.center = doCenter)
  
  
  seurat.obj <- RunPCA(object = seurat.obj,
                       features = VariableFeatures(object = seurat.obj),
                       verbose = FALSE, seed.use = 10, npc = npc)
  if(length(reg.var) > 0 ) seurat.obj = regress_on_pca(seurat.obj, reg.var)
  
  # seurat.obj <- RunLSI(seurat.obj, n = npc,
  #                      features = VariableFeatures(object = seurat.obj))
  
  return(seurat.obj)
}
doBasicSeurat_new = cmpfun(doBasicSeurat_new)


## query the resoltuion parameters given a seurat object and the number of clusters
## using binary seach
queryResolution4Seurat <- function(seurat.obj, k = 10, reduction = 'umap', npc = 20, 
                                   min_resl = 0.01, max_resl = 1, max_iter = 15, doPCA = F){
  max.dim = ifelse(reduction == 'pca', npc, 2)
  if(doPCA) {
    seurat.obj <- RunPCA(seurat.obj, npcs = npc, verbose = F)
    seurat.obj <- RunTSNE(seurat.obj, dims = 1:npc, verbose = F)
    seurat.obj <- RunUMAP(seurat.obj, dims = 1:npc, verbose = F)
  }
  
  
  seurat.obj <- FindNeighbors(seurat.obj, reduction = reduction, verbose = F, dims = 1:max.dim)
  tmp.cluster1 <- FindClusters(seurat.obj, resolution = min_resl)@active.ident
  tmp.cluster2 <- FindClusters(seurat.obj, resolution = max_resl)@active.ident
  
  
  
  
  len1 = length(levels(tmp.cluster1))
  len2 = length(levels(tmp.cluster2))
  
  k1 = k2 = 0
  while(len1 > k ){
    
    k1 = k1 + 1
    message('min_resl too large, trying to divided it by  2')
    min_resl = min_resl/2
    tmp.cluster1 <- FindClusters(seurat.obj, resolution = min_resl)@active.ident
    len1 = length(levels(tmp.cluster1))
    if(k1 == 10) return(min_resl)
    #if(k1 == 10) stop('Please specify a much smaller min_res')
  }
  
  while(len2 < k){
    k2 = k2 + 1
    message('max_resl too small, trying to multiply it by 2')
    max_resl = max_resl * 2
    tmp.cluster2 <- FindClusters(seurat.obj, resolution = max_resl)@active.ident
    len2 = length(levels(tmp.cluster2))
    #if(k2 == 10) stop('Please specify a much bigger max_res')
    if(k2 == 10) return(max_resl)
  }
  if(len1 == k) {
    return(min_resl)
  }
  
  if(len2 == k) {
    return(max_resl)
  }
  
  # repeat in other case
  
  i = 0
  repeat{
    i = i + 1
    resl0 = min_resl/2 + max_resl/2
    
    tmp.cluster <- FindClusters(seurat.obj, resolution = resl0)@active.ident
    
    len = length(levels(tmp.cluster)) 
    if(len == k){
      return(resl0)
    }
    if(len < k){
      min_resl = resl0
      len1 = len
    }
    if(len > k){
      max_resl = resl0
      len2 = len
    }
    if(i == max_iter) break
  }
  return(resl0)
}
queryResolution4Seurat = cmpfun(queryResolution4Seurat)

# do clustering using different #pcs and different methods, given an seurat.obj
## return seurat object, with a new meta.data column with name clusterLabelName
## gcSNN -- the clustering method used by seurat (Louvain algorithm on snn)
clust4seurat <- function(seurat.obj, npc = 30, method = 'gcSNN', reduction = 'tsne', resolution = 0.1, 
                         clustLabelName = paste0('clusterBy_', reduction, '_',  method), k = NULL){
  
  ## note the seurat object was done pca using 100 pcs; so npc should be smaller than 100
  ## using gcSNN
  if(is.null(seurat.obj@reductions$pca)){
    seurat.obj <- RunPCA(seurat.obj, npcs = npc, verbose = F)
  } 
  
  
  if(reduction == 'tsne'){
    seurat.obj <- RunTSNE(seurat.obj, dims = 1:npc)
  }
  if(reduction == 'umap'){
    seurat.obj <- RunUMAP(seurat.obj, dims = 1:npc, verbose = F)
  }
  
  
  
  if(method == 'gcSNN'){
    
    if(!is.null(k)){
      # find best resolution to get k cluster
      resolution = queryResolution4Seurat(seurat.obj, reduction = reduction, npc = npc, min_resl = resolution,
                                          max_resl = 10*resolution, k = k)
    }
    
    if(reduction == 'tsne'){
      
      seurat.obj <- FindNeighbors(seurat.obj, reduction = reduction, dims = 1:2)
      seurat.obj <- FindClusters(seurat.obj, resolution = resolution, verbose = F)
    }
    
    if(reduction == 'umap'){
      seurat.obj <- FindNeighbors(seurat.obj, reduction = reduction, dims = 1:2)
      seurat.obj <- FindClusters(seurat.obj, resolution = resolution, verbose = F)
    }
    
    if(reduction == 'pca'){
      seurat.obj <- FindNeighbors(seurat.obj, reduction = reduction, dims = 1:npc)
      seurat.obj <- FindClusters(seurat.obj, resolution = resolution, verbose = F)
    }
    seurat.obj@meta.data[[clustLabelName]] = as.integer(seurat.obj@active.ident)
    
  }else{
    if(reduction == 'tsne'){
      
      embd.mat <- seurat.obj@reductions$tsne@cell.embeddings
      
    }
    
    if(reduction == 'umap'){
      
      embd.mat <- seurat.obj@reductions$umap@cell.embeddings
    }
    
    if(reduction == 'pca'){
      embd.mat <- seurat.obj@reductions$pca@cell.embeddings[, 1:npc]
    }
    
    
    seurat.obj@meta.data[[clustLabelName]] = generalCluster(embd.mat, method, k)
  }
  
  
  
  return(seurat.obj)
}



## do graph cluster snn on cistop object
gcSNN4CistopicObj <- function(cisTopic.obj, ntopic = 20, clust_by = 'topic',
                              resolution = 0.1, k = NULL){
  #note cisTopic.obj should include model with ntopic
  sele.cisTopic <- selectModel(cisTopic.obj, select = ntopic, 
                               keepBinaryMatrix = F, keepModels = F)
  sele.cisTopic <- runtSNE(sele.cisTopic, target='cell')
  sele.cisTopic <- runUmap(sele.cisTopic, target='cell')
  
  if(clust_by == 'topic'){
    embd.mtx = t(modelMatSelection(sele.cisTopic, 'cell', 'Probability'))
  }
  if(clust_by == 'tsne'){
    embd.mtx = sele.cisTopic@dr$cell$tSNE
  }
  if(clust_by == 'umap'){
    embd.mtx = sele.cisTopic@dr$cell$Umap
  }
  
  snn = FindNeighbors(dist(embd.mtx))
  names(snn) = c('nn', 'snn')
  # creat a seurat object to faciliat finding clusters using seurat
  tmp.seurat <- CreateSeuratObject(sele.cisTopic@count.matrix)
  tmp.seurat@graphs = snn
  if(!is.null(k)) {
    
    resolution = queryResolution4Topic(tmp.seurat, k = k, min_resl = resolution, max_resl = 10 * resolution)
  }
  
  tmp.seurat <- FindClusters(tmp.seurat, graph.name = 'snn', resolution = resolution)
  #tmp.seurat@reductions$tsne@cell.embeddings = sele.cisTopic@dr$cell$tSNE
  #tmp.seurat@reductions$umap@cell.embeddings = sele.cisTopic@dr$cell$Umap
  
  colnames(sele.cisTopic@dr$cell$tSNE) = c('tSNE_1', 'tSNE_2')
  colnames(sele.cisTopic@dr$cell$Umap) = c('UMAP_1', 'UMAP_2')
  
  return(list('tsne_coord' = sele.cisTopic@dr$cell$tSNE, 
              'umap_coord' = sele.cisTopic@dr$cell$Umap,
              'cluster_label' = as.integer(tmp.seurat@active.ident)))
}


## run sc3
run_sc3 <- function(mtx, k = 5){
  fnames = rownames(mtx)
  rownames(mtx) = NULL
  sce <- SingleCellExperiment(
    assays = list(
      counts = as.matrix(mtx),
      logcounts = log2(as.matrix(mtx) + 1)
    )
  )
  rowData(sce)$feature_symbol = fnames
  
  sc3 = sc3(sce, ks = k, biology = T)
  
  return(colData(sc3))
}


## run monocle3
run_monocle3 <- function(mtx, nFeatures = 10000){
  pd = data.frame('cellID' = 1:ncol(mtx))
  row.names(pd) = colnames(mtx)
  
  
  fd = data.frame('gene_short_name' = rownames(mtx))
  rownames(fd) = rownames(mtx)
  
  pd <- new("AnnotatedDataFrame", data = pd)
  fd <- new("AnnotatedDataFrame", data = fd)
  
  mtx.monocle <- newCellDataSet(mtx,
                                phenoData = pd,
                                featureData = fd,
                                expressionFamily=VGAM::negbinomial.size(),
                                lowerDetectionLimit=1)
  
  mtx.monocle <- estimateSizeFactors(mtx.monocle)
  mtx.monocle <- estimateDispersions(mtx.monocle)
  
  
  
  disp_table = dispersionTable(mtx.monocle)
  disp_table = disp_table %>% mutate(excess_disp =
                                       (dispersion_empirical - dispersion_fit) / dispersion_fit) %>%
    arrange(plyr::desc(excess_disp))
  
  top_subset_genes = as.character(head(disp_table, nFeatures)$gene_id)
  
  mtx.monocle = setOrderingFilter(mtx.monocle, top_subset_genes)
  mtx.monocle <- preprocessCDS(mtx.monocle,  method = 'PCA',
                               norm_method = 'log',
                               num_dim = 50,
                               verbose = F)
  
  if(F){
    mtx.monocle <- reduceDimension(mtx.monocle, max_components = 2,
                                   reduction_method = 'UMAP',
                                   metric="correlation",
                                   min_dist = 0.75,
                                   n_neighbors = 50,
                                   verbose = T)
    
    mtx.monocle <- clusterCells(mtx.monocle,
                                method = 'louvain',
                                louvain_iter = 1,
                                verbose = T)
    
  }
  
  #plot_cell_clusters(mtx.monocle,
  #                   color_by = 'Cluster',
  #                   cell_size = 0.1,
  #                   show_group_id = T)  
  
  
  #meta.data = pData(mtx.monocle)
  return(mtx.monocle)
}


## implement of scATAC-seq sepecific methods ####
# the imput mtx is already filterd
run_scABC <- function(mtx, k = 5){
  weights = apply(mtx, 2, mean)
  landmarks = computeLandmarks(mtx, weights = weights, nCluster = k)
  labels = assign2landmarks(mtx, landmarks)
  return(labels)
}


run_chromVAR <- function(mtx, genomeName = 'BSgenome.Hsapiens.UCSC.hg19'){
  
  register(MulticoreParam(4))
  rs = Matrix::rowSums(mtx)
  mtx = mtx[rs>0, ] 
  peaks = data.table('x' = rownames(mtx))
  peaks = tidyr::separate(peaks, col = 'x', into = c('chr', 'start', 'end'))
  peaks = GenomicRanges::makeGRangesFromDataFrame(peaks)
  
  frag.counts = SummarizedExperiment(assay = list(counts = mtx),
                                     rowRanges = peaks)
  frag.counts <- addGCBias(frag.counts, genome = genomeName)
  motifs <- getJasparMotifs()
  motif_ix <- matchMotifs(motifs, frag.counts,
                          genome = genomeName)
  dev <- computeDeviations(object = frag.counts, 
                           annotations = motif_ix)
  bg <- getBackgroundPeaks(object = frag.counts)
  
  dev <- computeDeviations(object = frag.counts, annotations = motif_ix,
                           background_peaks = bg)
  
  #motif.zscore = dev@assays$data$z
  return(dev)
}


run_LSI <- function(mtx, ncell.peak = 150,  max_pc = 10, k = 5){
  mtx = (mtx > 0)
  mtx = 1 * mtx
  
  num_cells.peak = Matrix::rowSums(mtx)
  ncounts = mtx[num_cells.peak >= ncell.peak,]
  
  ## Normalize the data with TF-IDF
  nfreqs = t(t(ncounts) / Matrix::colSums(ncounts))
  tf_idf_counts = nfreqs * log(1 + ncol(ncounts) / Matrix::rowSums(ncounts))
                               
  ## DO SVD on tf_idf normalized matrix                             
  set.seed(0)
  SVDtsne = irlba(tf_idf_counts, max_pc, max_pc, maxit=1000)
  d_diagtsne = matrix(0, nrow=length(SVDtsne$d), ncol=length(SVDtsne$d))
  diag(d_diagtsne) = SVDtsne$d
  SVDtsne_vd = t(d_diagtsne %*% t(SVDtsne$v))
  rownames(SVDtsne_vd) = colnames(mtx)
  colnames(SVDtsne_vd) = paste0('pca_', 1:ncol(SVDtsne_vd))
  
  ## Run TSNE to 2 dimensions
  if(F){
    tsnetfidf = Rtsne(SVDtsne_vd, pca=F, perplexity=30, max_iter=5000)
    
    tsne_coords = as.data.frame(tsnetfidf$Y)
    colnames(tsne_coords) = c('tsne_1', 'tsne_2')
    rownames(tsne_coords) = colnames(ncounts) 
  }
   
  pca_coords = SVDtsne_vd[, 2:max_pc]
  cl.labels = generalCluster(pca_coords, k = k, method = 'hclust')
  return(cl.labels)
}


run_scrat <- function(mtx, reduction = 'pca', max_pc = 20, method = 'mclust', k = 10){
  # stadardized features per cell
  mtx = scale(mtx, center = T, scale = T)
  if(reduction == 'pca'){
    reduced.mtx = doDimReduction4mat(mtx, max_pc = max_pc)[[1]]
    reduced.mtx = reduced.mtx[, -1]
  }
  if(reduction == 'tsne'){
    reduced.mtx = doDimReduction4mat(mtx, max_pc = max_pc, doTSNE = T)[[2]]
  }
  
  cl.label = generalCluster(reduced.mtx, k = k, method = method)
  return(cl.label)
}

# fit cistopic model and using the predicted cell * topic probability matrix
# for clustering
run_cisTopic <- function(mtx, nCores = 4){
  # prepare the right format of rownames
  rnames = data.table('region' = rownames(mtx))
  tmp = tidyr::separate(rnames, col = 'region', into = c('chr', 'start', 'end'))
  rnames = paste0(tmp$chr, ':', tmp$start, '-', tmp$end)
  rownames(mtx) = rnames
  
  cisTopicObject <- createcisTopicObject(mtx, project.name='scATAC')
  cisTopicObject <- runModels(cisTopicObject, topic = c(5, 10, 15, 20, 30, 40,
                                                        50, 60, 70, 80, 90, 100), seed = 987, nCores = nCores, 
                              burnin = 120, iterations = 150, addModels = T)
  #cisTopicObject <- selectModel(cisTopicObject, keepBinarymatrix = F, keepModels = F)
  #cellassign <- t(modelMatSelection(cisTopicObject, 'cell', 'Probability'))
  return(cisTopicObject)
}



queryResolution4Monocle <- function(obj, k = 5, 
                                    min_resl = 0.1, max_resl = 1, knn = 50,
                                    max_iter = 15, doPCA = F){
  
  tmp.cluster1 <- clusterCells(obj, method = 'louvain',
                      louvain_iter = 1, res = min_resl * 10^(-2),
                      verbose = F, k = knn)
  tmp.cluster2 <- clusterCells(obj, method = 'louvain',
                       louvain_iter = 1, res = max_resl * 10^(-2),
                       verbose = F, k = knn)

  len1 = length(levels(pData(tmp.cluster1)$Cluster))
  len2 = length(levels(pData(tmp.cluster2)$Cluster))
  
  k1 = k2 = 0
  while(len1 > k ){
    
    k1 = k1 + 1
    message('min_resl too large, trying to divided it by  2')
    min_resl = min_resl/2
    tmp.cluster1 <- clusterCells(obj, method = 'louvain',
                                 louvain_iter = 1, res = min_resl * 10^(-2),
                                 verbose = F, k = knn)
    len1 = length(levels(pData(tmp.cluster1)$Cluster))
    if(k1 == 10) stop('Please specify a much smaller min_res')
  }
  
  while(len2 < k){
    k2 = k2 + 1
    message('max_resl too small, trying to multiply it by 2')
    max_resl = max_resl * 2
    tmp.cluster2 <- clusterCells(obj, method = 'louvain',
                                 louvain_iter = 1, res = max_resl * 10^(-2),
                                 verbose = F, k = knn)
    len2 = length(levels(pData(tmp.cluster2)$Cluster))
    if(k2 == 10) stop('Please specify a much bigger max_res')
  }
  if(len1 == k) {
    return(min_resl)
  }
  
  if(len2 == k) {
    return(max_resl)
  }
  
  # repeat in other case
  
  i = 0
  repeat{
    i = i + 1
    resl0 = min_resl/2 + max_resl/2
    
    tmp.cluster <- clusterCells(obj, method = 'louvain',
                                louvain_iter = 1, res = resl0 * 10^(-2),
                                verbose = F, k = knn)
    
    len = length(levels(pData(tmp.cluster)$Cluster))
    if(len == k){
      return(resl0)
    }
    if(len < k){
      min_resl = resl0
      len1 = len
    }
    if(len > k){
      max_resl = resl0
      len2 = len
    }
    if(i == max_iter) break
  }
  return(resl0)
}



## implement of generic clustering methods ####
generalCluster <- function(reduced.mtx, method = 'hclust', k = 5){
  if(method == 'kmeans'){
    res = kmeans(reduced.mtx, centers = k)
    cl.label = res$cluster
  }
  
  if(method == 'hclust'){
    if(is.null(k)) stop('Need specify k: the number of cluster')
    d <- dist(reduced.mtx, method = "euclidean") # distance matrix
    fit <- hclust(d, method = "ward.D")
    cl.label<- cutree(fit, k = k) # cut tree into k clusters
  }
  
  if(method == 'mclust'){
    
    fit <- Mclust(reduced.mtx, G = k)
    cl.label = fit$classification
  }
  return(cl.label)
}



