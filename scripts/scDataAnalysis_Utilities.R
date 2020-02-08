
## source functions used for analysis singel-cell data in unit of TAD

## list of functions ####
library(data.table)
library(compiler)
library(magrittr)
library(Rtsne)
library(Seurat)
library(DropletUtils)
library(cicero)
library(cisTopic)
library(ggpubr)
library(RColorBrewer)
#library(scABC)
#library(preprocessCore)

# for clustering and/or validation
library(cluster)
library(factoextra)
library(fpc)
library(mclust)

## **** note: for cluster validation, check:  ***** ##
## **** http://www.sthda.com/english/wiki/wiki.php?id_contents=7952  ***##


## do reverse complemente of a DNA sequence
rev.comp <- function(x, rev=TRUE){
  x<-toupper(x)
  y<-rep("N",nchar(x))
  xx<-unlist(strsplit(x,NULL))
  for (bbb in 1:nchar(x))
  {
    if(xx[bbb]=="A") y[bbb]<-"T"    
    if(xx[bbb]=="C") y[bbb]<-"G"    
    if(xx[bbb]=="G") y[bbb]<-"C"    
    if(xx[bbb]=="T") y[bbb]<-"A"
  }
  if(rev==FALSE) 
  {
    for(ccc in (1:nchar(x)))
    {
      if(ccc==1) yy<-y[ccc] else yy<-paste(yy,y[ccc],sep="")
    }
  }
  if(rev==T)
  {
    zz<-rep(NA,nchar(x))
    for(ccc in (1:nchar(x)))
    {
      zz[ccc]<-y[nchar(x)+1-ccc]
      if(ccc==1) yy<-zz[ccc] else yy<-paste(yy,zz[ccc],sep="")
    }
  }
  return(yy)  
}
# rebin data in unit of tad, adjusting tad size

rebin_matrix2tad <- function(mtx, tads){
  # mtx: matrix wiht rownames as chr-start-end and
  # colnames as cell names
  tads[, 'size' := floor(end/1000 - start/1000)]
  setkey(tads, id)
  
  rnames = rownames(mtx)
  mtx_chr = sapply(rnames, function(x) unlist(strsplit(x, '-'))[1])
  chrs = unique(mtx_chr)
  starts = as.numeric(sapply(rnames, function(x) unlist(strsplit(x, '-'))[2]))
  ends = as.numeric(sapply(rnames, function(x) unlist(strsplit(x, '-'))[3]))
  
  rebin_mat = NULL
  for(chr0 in chrs){
    mtx0 = mtx[mtx_chr == chr0, ]
    mtx0 = as.data.table(mtx0)
    mtx0$start = starts[mtx_chr == chr0] 
    mtx0$end = ends[mtx_chr == chr0] 
    mtx0[, 'midP' := start/2 + end/2]
    tads0 = tads[chr == chr0]
    mtx0[, 'tad_id' := which(tads0$start <= midP & tads0$end >= midP), by = midP]
    mtx0[, tad_id := tads0$id[tad_id]]
    mtx0[, c('start', 'end', 'midP') := NULL]
    rebin_mat = rbind(rebin_mat, mtx0)
  }
  
  rebin_mat = data.table(rebin_mat)
  setkey(rebin_mat, tad_id)
  
  new_mat = rebin_mat[, lapply(.SD, sum), by = tad_id]
  new_mat = new_mat[complete.cases(new_mat)]
  
  feature.names = new_mat$tad_id
  new_mat[, 'size' := tads[J(new_mat$tad_id)]$size]
  
  
  new_mat = new_mat[, lapply(.SD, function(x) x/size * 1000),
                    .SDcols = !c('size', 'tad_id')]
  
  new_mat = as.matrix(new_mat)
  rownames(new_mat) = feature.names
  
  return(new_mat)
}
rebin_matrix2tad = cmpfun(rebin_matrix2tad)


rebin_matrix2Bin <- function(mtx, resl = 100 * 1000){
  # mtx: matrix wiht rownames as chr-start-end and
  # colnames as cell names
  
  rnames = rownames(mtx)
  mtx_chr = sapply(rnames, function(x) unlist(strsplit(x, '-'))[1])
  chrs = unique(mtx_chr)
  starts = as.numeric(sapply(rnames, function(x) unlist(strsplit(x, '-'))[2]))
  ends = as.numeric(sapply(rnames, function(x) unlist(strsplit(x, '-'))[3]))
  
  rebin_mat = NULL
  for(chr0 in chrs){
    mtx0 = mtx[mtx_chr == chr0, ]
    mtx0 = as.data.table(mtx0)
    mtx0$start = starts[mtx_chr == chr0] 
    mtx0$end = ends[mtx_chr == chr0] 
    
    mtx0[, 'id' := ceiling((start+ end)/resl/2)]
    mtx0[, 'bin_id' := paste0(chr0, '-', id)]
    mtx0[, c('start', 'end', 'id') := NULL]
    rebin_mat = rbind(rebin_mat, mtx0)
  }
  
  rebin_mat = data.table(rebin_mat)
  setkey(rebin_mat, bin_id)
  
  new_mat = rebin_mat[, lapply(.SD, sum), by = bin_id]
  new_mat = new_mat[complete.cases(new_mat)]
  
  feature.names = new_mat$bin_id
  new_mat[, 'bin_id' := NULL]
  
  
  new_mat = as.matrix(new_mat)
  rownames(new_mat) = feature.names
  
  return(new_mat)
}
rebin_matrix2Bin = cmpfun(rebin_matrix2Bin)


segByPeak <- function(peaks){
  peaks[, 'midP' := floor(start/2 + end/2)]
  chrs = unique(peaks$chr)
  domains = NULL
  for(chr0 in chrs){
    bd0 = peaks[chr == chr0]
    setkey(bd0, midP)
    len = nrow(bd0)
    tmp = data.table('chr' = chr0, 'start' = bd0$midP[1:(len - 1)] - 500, 
                     'end' = bd0$midP[2:len] + 500)
    domains = rbind(domains, tmp)
  }
  return(domains)
}
segByPeak = cmpfun(segByPeak)

## remove gap between tads
full_seg_tads <- function(tads){
  chrs = unique(tads$chr)
  res = NULL
  for(chr0 in chrs){
    tads0 = tads[chr == chr0]
    bounds = sort(unique(c(tads0$start, tads0$end)))
    len = length(bounds)
    res = rbind(res, data.table('chr' = chr0,
                                'start' = bounds[-len], 'end' = bounds[-1]))
  }
  
  return(res)
}
full_seg_tads = cmpfun(full_seg_tads)


# evaluate clustering given cell*feature matrix object
eval_cluster4mat <- function(mat, cluster.label, alt.cluster.label = NULL, distMethod = 'euclidean'){
 
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

# evaluate clustering given seurat object
eval_cluster4seurat <- function(seurat.obj, reduction = 'tsne', npc = 30,
                                distMethod = 'euclidean', alt.cluster.label = NULL){
  if(reduction == 'tsne'){
    label1 =  seurat.obj@active.ident
    mat = seurat.obj@reductions$tsne@cell.embeddings
  }
  
  if(reduction == 'umap'){
    label1 =  seurat.obj@active.ident
    mat = seurat.obj@reductions$umap@cell.embeddings
  }
  
  if(reduction == 'pca'){
    label1 =  seurat.obj@active.ident
    mat = seurat.obj@reductions$pca@cell.embeddings[, 1:npc]
  }
  
  Dist = dist(mat, method = distMethod)
  
  label1 = as.numeric(label1)
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

# do normalization, pca using Seurat
doBasicSeurate <- function(mtx, npc = 50, top.variable = 0.5, doLog = T, 
                           doScale = T, doCenter = T, assay = 'ATAC', 
                           regr.var = 'nCount_ATAC'){
  
 # top.variabl -- use top most variable features
  if(doLog) mtx = log2(1 + mtx)
  seurat.obj = CreateSeuratObject(mtx, project = 'scATAC', assay = assay,
                                  names.delim = '-')
  cell.names = colnames(mtx)
  
  #seurat.obj <- NormalizeData(seurat.obj, normalization.method = 'LogNormalize',
  #                            scale.factor = 1e4)
  
  seurat.obj <- FindVariableFeatures(object = seurat.obj, 
                                     selection.method = 'vst', 
                                     nfeatures = floor(nrow(mtx) * top.variable))
  seurat.obj <- ScaleData(object = seurat.obj, 
                          features = VariableFeatures(seurat.obj), 
                          vars.to.regress = regr.var, do.scale = doScale, do.center = doCenter)
  
  
  seurat.obj <- RunPCA(object = seurat.obj, 
                       features = VariableFeatures(object = seurat.obj),
                       verbose = FALSE, seed.use = 10, npc = npc)
  
  return(seurat.obj)
}
doBasicSeurate = cmpfun(doBasicSeurate)


# do normalization, pca using Seurat
doBasicSeurate_RNA <- function(mtx, npc = 50, top.variable = 0.2, pmito.upper = 0.2,
                           doScale = T, doCenter = T, assay = 'RNA', 
                           regr.var = 'nCount_RNA'){

 ## top.variabl -- use top most variable features

 # filter cells with high percentage of mitocondria genes

  mito.features <- grep(pattern = "^MT-", 
                      x = rownames(x = seurat.obj), value = TRUE)

  perc.mito = Matrix::colSums(mtx[mito.features, ])/Matrix::colSums(mtx)

  mtx = mtx[, perc.mito <= pmito.upper]
  perc.mito = perc.mito[perc.mito <= pmito.upper]


 # create seurat object
  seurat.obj = CreateSeuratObject(mtx, project = 'scRNA', assay = assay,
                                  names.delim = '-', min.cells = 10, min.features = 100)
  

  if(all(names(perc.mito) == conames(mtx))) seurat.obj@meta.data[['perc.mito']] = perc.mito

  seurat.obj <- subset(x = seurat.obj, subset = (nFeature_RNA < 10000))

  
  seurat.obj <- NormalizeData(seurat.obj, normalization.method = 'LogNormalize',
                              scale.factor = 1e4)
  
  seurat.obj <- FindVariableFeatures(object = seurat.obj, 
                                     selection.method = 'vst', 
                                     nfeatures = floor(nrow(mtx) * top.variable))
  seurat.obj <- ScaleData(object = seurat.obj, 
                          features = VariableFeatures(seurat.obj), 
                          vars.to.regress = regr.var, do.scale = doScale, do.center = doCenter)
  
  
  seurat.obj <- RunPCA(object = seurat.obj, 
                       features = VariableFeatures(object = seurat.obj),
                       verbose = FALSE, seed.use = 10, npc = npc)
  
  return(seurat.obj)
}
doBasicSeurate_RNA = cmpfun(doBasicSeurate_RNA)



# in tad
doBasicSeurate_tad <- function(tads, mtx, npc = 50){
  
  tads = full_seg_tads(tads)
  tads[, 'id' := paste(chr, start, end, sep = '-')]
  rebinned_mtx = rebin_matrix2tad(mtx, tads)
  rebinned_mtx = log2(1 + rebinned_mtx)
  seurat.obj = CreateSeuratObject(rebinned_mtx, project = 'scATAC_tad', assay = 'ATAC',
                                  names.delim = '-')
  cell.names = colnames(rebinned_mtx)
  
  #seurat.obj <- NormalizeData(seurat.obj, normalization.method = 'LogNormalize',
  #                            scale.factor = 1e4)
 
  seurat.obj <- FindVariableFeatures(object = seurat.obj, 
                                     selection.method = 'dispersion', 
                                     nfeatures = nrow(rebinned_mtx))
  seurat.obj <- ScaleData(object = seurat.obj, 
                          features = VariableFeatures(seurat.obj), 
                          vars.to.regress = c('nCount_ATAC'), do.scale = T, do.center = T)
  
  
  seurat.obj <- RunPCA(object = seurat.obj, 
                       features = VariableFeatures(object = seurat.obj),
                       verbose = FALSE, seed.use = 10, npc = npc)
  
  return(seurat.obj)
}
doBasicSeurate_tad = cmpfun(doBasicSeurate_tad)

# in bin
doBasicSeurate_bin <- function(mtx, resl = 500 * 1000, npc = 50){
  
  rebinned_mtx = rebin_matrix2Bin(mtx, resl)
  rebinned_mtx = log2(1 + rebinned_mtx)
  seurat.obj = CreateSeuratObject(rebinned_mtx, project = 'scATAC_tad', assay = 'ATAC',
                                  names.delim = '-')
  cell.names = colnames(rebinned_mtx)
  
  seurat.obj <- FindVariableFeatures(object = seurat.obj, 
                                     selection.method = 'dispersion', 
                                     nfeatures = nrow(rebinned_mtx))
  seurat.obj <- ScaleData(object = seurat.obj, 
                          features = VariableFeatures(seurat.obj), 
                          vars.to.regress = c('nCount_ATAC'), do.scale = T, do.center = T)
  
  
  seurat.obj <- RunPCA(object = seurat.obj, 
                       features = VariableFeatures(object = seurat.obj),
                       verbose = FALSE, seed.use = 10, npc = npc)
  
  return(seurat.obj)
}
doBasicSeurate_bin = cmpfun(doBasicSeurate_bin)


# in given ctcf peaks
doBasicSeurate_ctcf <- function(ctcf_peaks, mtx, npc = 50){
  
  ptads = segByPeak(ctcf_peaks)
  ptads[, 'id' := paste(chr, start, end, sep = '-')]
  ptads[, 'size' := floor((end - start)/1000)]
  ptads = ptads[size > 2 & size < 3000]
  rebinned_mtx = rebin_matrix2tad(mtx, ptads)
  rebinned_mtx = log2(1 + rebinned_mtx)
  seurat.obj = CreateSeuratObject(rebinned_mtx, project = 'scATAC_tad', assay = 'ATAC',
                                  names.delim = '-')
  cell.names = colnames(rebinned_mtx)
  
  #seurat.obj <- NormalizeData(seurat.obj, normalization.method = 'LogNormalize',
  #                            scale.factor = 1e4)
  
  seurat.obj <- FindVariableFeatures(object = seurat.obj, 
                                     selection.method = 'dispersion', 
                                     nfeatures = nrow(rebinned_mtx) * 0.5)
  seurat.obj <- ScaleData(object = seurat.obj, 
                          features = VariableFeatures(seurat.obj), 
                          vars.to.regress = c('nCount_ATAC'), do.scale = T, do.center = T)
  
  
  seurat.obj <- RunPCA(object = seurat.obj, 
                       features = VariableFeatures(object = seurat.obj),
                       verbose = FALSE, seed.use = 10, npc = npc)
  
  return(seurat.obj)
}
doBasicSeurate_ctcf = cmpfun(doBasicSeurate_ctcf)


# basic seurat plot
basicSeuratPlot <- function(org.seurat.obj){
  p1 <- VizDimLoadings(object = org.seurat.obj, dims = 1:2, nfeatures = 20)
  
  p2 <- DimHeatmap(object = org.seurat.obj, dims = 1:6, cells = 100,   balanced = TRUE)
  p3 <- ElbowPlot(object = org.seurat.obj)
  return(list(p1, p2, p3))
}

# cluster on tsne/umap/pca
doClusterBySeurate <- function(seurat.obj, reduction = 'tsne', npc = 30,
                               dims = 1:2, resolution = 0.6, start_pc = 1){
  if(reduction == 'tsne'){
    seurat.obj <- RunTSNE(object = seurat.obj, dims = start_pc:npc, reduction = 'pca')
    seurat.obj <- FindNeighbors(object = seurat.obj, dims = dims, reduction = 'tsne',
                                verbose = F)
  }
  if(reduction == 'umap'){
    seurat.obj <- RunUMAP(object = seurat.obj, dims = start_pc:npc, reduction = 'pca', verbose = F)
    seurat.obj <- FindNeighbors(object = seurat.obj, dims = dims, reduction = 'umap',
                                verbose = F)
  }
  if(reduction == 'pca'){

    seurat.obj <- FindNeighbors(object = seurat.obj, dims = dims, reduction = 'pca',
                                verbose = F)
    seurat.obj <- RunTSNE(object = seurat.obj, dims = start_pc:npc, reduction = 'pca')
    seurat.obj <- RunUMAP(object = seurat.obj, dims = start_pc:npc, reduction = 'pca', verbose = F)
  }
  
  seurat.obj <- FindClusters(object = seurat.obj, resolution = resolution, verbose = F)
  
  return(seurat.obj)
}
doClusterBySeurate = cmpfun(doClusterBySeurate)


# integrate analysis
# data were given in unit of TAD
doSeurat_integrate_tad <- function(atac.mtx, gene.mtx, qLorm = T){
  colnames(atac.mtx) = paste0('atac-', colnames(atac.mtx))
  comb.mtx = cbind(gene.mtx, atac.mtx)
  dtype = sapply(colnames(comb.mtx), function(x) ifelse(grepl(x, pattern = '^atac-'), 'ATAC', 'RNA'))
  
  if(qLorm) {
    rnames = rownames(comb.mtx)
    cnames = colnames(comb.mtx)
    rand.mat = matrix(runif(nrow(comb.mtx) * ncol(comb.mtx), 0, 10^(-6)), nrow(comb.mtx))
    comb.mtx <- normalize.quantiles(as.matrix(comb.mtx) + rand.mat)
    colnames(comb.mtx) = cnames
    rownames(comb.mtx) = rnames
  }
  seurat.integrate <- CreateSeuratObject(comb.mtx, assay = 'Comb')
  seurat.integrate@meta.data$orig.ident = dtype

  #seurat.obj <- NormalizeData(seurat.integrate, normalization.method = 'LogNormalize',
  #                            scale.factor = 1e4)
  
  seurat.integrate <- FindVariableFeatures(object = seurat.integrate, 
                                     selection.method = 'dispersion', 
                                     nfeatures = nrow(comb.mtx))
  seurat.integrate <- ScaleData(object = seurat.integrate, 
                          features = VariableFeatures(seurat.integrate), 
                          vars.to.regress = c('orig.ident', 'nCount_Comb'), do.scale = F, do.center = F)
  
  
  seurat.integrate <- RunPCA(object = seurat.integrate, 
                       features = VariableFeatures(object = seurat.integrate),
                       verbose = FALSE, seed.use = 10, npc = 50)
  return(seurat.integrate)
}


## map gene to atac peak
gene2peak <- function(gene_set, peaks, gene_ann){
  # should include tss information in gene_list
  gene_list = gene_ann[gene_name %in% gene_set, ]
  chrs = unique(gene_list$chr)
  gene_new = NULL
  peaks[, 'midP' := start/2 + end/2]
  for(chr0 in chrs){
    gene0 = gene_list[chr == chr0, ]
    peaks0 = peaks[chr == chr0]
    gene0[, 'peak_id' := which(tss >= peaks0$start & tss <= peaks0$end), by = gene_id]
    gene0[, 'peak_id' := ifelse(is.na(peak_id), which.min(abs(tss - peaks0$midP)), peak_id), by = gene_name]
    gene0[, 'peak' := peaks0[peak_id]$pos]
    
    gene_new = rbind(gene_new, gene0)
  }
  return(gene_new)
}


# plot atac signal around Tss, given a gene set and seurat atac object
plotSignal_Tss <- function(genes, org.seurat.obj, reduction = 'tsne'){
  ## genes should include gene_name and the corresponding peak information
  tmp_plot = list()
  for(i in 1:nrow(genes)){
    tmp_plot[[i]] <- FeaturePlot(object = org.seurat.obj, reduction = reduction,  feature = genes$peak[i]) + labs(title = genes$gene_name[i])
  } 
  
  pp <- ggarrange(plotlist = tmp_plot, nrow = 2, ncol = ceiling(length(tmp_plot)/2))
  return(pp)
}

read10X_ATAC <- function(dirt){
  mtx_path <- paste0(dirt, "matrix.mtx")
  feature_path <- paste0(dirt, "peaks.bed")
  barcode_path <- paste0(dirt, "barcodes.tsv")
  
  
  features <-readr::read_tsv(feature_path, col_names = F) %>% tidyr::unite(feature, sep = '-')
  barcodes <- readr::read_tsv(barcode_path, col_names = F) %>% tidyr::unite(barcode)
  
  mtx <-  Matrix::readMM(mtx_path) %>%
    magrittr::set_rownames(features$feature)%>%
    magrittr::set_colnames(barcodes$barcode) 
  
  return(mtx)
}

compare_cellCall4ATAC <- function(cellRanger_dir, fdr = 0.001, lower = 500, upper = NULL){
  filter_dir = paste0(cellRanger_dir, 'filtered_peak_bc_matrix/')
  raw_dir = paste0(cellRanger_dir, 'raw_peak_bc_matrix/')
  
  mat <- read10X_ATAC(raw_dir)
  
  
  set.seed(2019)
  cell.out <- emptyDrops(mat, lower = lower, retain = upper)
  
  filter.out <- cell.out[complete.cases(cell.out), ]
  
  #is.cell <- (cell.out$FDR <= fdr)
  #plot(cell.out$Total, -cell.out$LogProb, col=ifelse(is.cell, "red", "black"),
  #     xlab="Total count", ylab="-Log Probability", log = 'x')
  
  
  filter.out = filter.out[filter.out$FDR <= fdr, ]
  
  rm(mat)
  
  
  mat <- read10X_ATAC(filter_dir)
  
  overlapped.cell <- intersect(rownames(filter.out), colnames(mat))
  dim(mat)
  perc.in.cellranger <- length(overlapped.cell)/nrow(filter.out)
  perc.in.emptydrop <- length(overlapped.cell)/ncol(mat)
  
  
  
  output = list('emptyDrop_cells' = filter.out, 'perc.in.cellranger' = round(perc.in.cellranger, 3),
                'perc.in.emptydrop' = round(perc.in.emptydrop, 3))
  
  return(output)
}


# do cicero given a seurate object
doCicero <- function(seurat.obj, reduction = 'tsne', chr_sizes,
                     gene_ann, npc = 30){
  ## gene_ann: the first four columns: chr, start, end, gene name
  set.seed(2019) 

  mtx = GetAssayData(seurat.obj, slot = 'counts')
  # change rownames using _ to delimited
  rnames = rownames(mtx)
  new.rnames = sapply(rnames, function(x) gsub('-', '_', x))
  rownames(mtx) <- new.rnames
  
  dt = reshape2::melt(as.matrix(mtx), value.name = 'count')
  rm(mtx)
  dt = dt[dt$count > 0, ]
  input_cds <- make_atac_cds(dt, binarize = T)
  rm(dt)
  input_cds <- detectGenes(input_cds)
  input_cds <- estimateSizeFactors(input_cds)
  
  if(reduction == 'tsne') {
    if(is.null(seurat.obj@reductions$tsne))
      seurat.obj <- RunTSNE(seurat.obj, dims = 1:npc)
    redu.coords = seurat.obj@reductions$tsne@cell.embeddings
  }
  if(reduction == 'umap') {
    if(is.null(seurat.obj@reductions$umap))
      seurat.object <- RunUMAP(seurat.object, dims = 1:npc)
    redu.coords = seurat.obj@reductions$umap@cell.embeddings
  }
  
  #make the cell id consistet
  
  cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = redu.coords)
  
  ## get connections
  
  conns <- run_cicero(cicero_cds, chr_sizes) 
  
  ## get cicero gene activity score
  names(gene_ann)[4] <- "gene"
  
  input_cds <- annotate_cds_by_site(input_cds, gene_ann)
  
  # generate unnormalized gene activity matrix
  unnorm_ga <- build_gene_activity_matrix(input_cds, conns)
  
  # make a list of num_genes_expressed
  num_genes <- pData(input_cds)$num_genes_expressed
  names(num_genes) <- row.names(pData(input_cds))
  
  # normalize
  cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)
  
  # if you had two datasets to normalize, you would pass both:
  # num_genes should then include all cells from both sets
  #unnorm_ga2 <- unnorm_ga
  #cicero_gene_activities <- normalize_gene_activities(list(unnorm_ga, unnorm_ga2), num_genes)
  
  #res = list('conns' = conns, 'gene_activities' = cicero_gene_activities)
  return(cicero_gene_activities)
}


# do cicero given a seurate object, just return the connection 
doCicero_modify <- function(seurat.obj, reduction = 'tsne', chr_sizes, npc = 30){
  ## gene_ann: the first four columns: chr, start, end, gene name
  set.seed(2019) 

  mtx = GetAssayData(seurat.obj, slot = 'counts')
  # change rownames using _ to delimited
  rnames = rownames(mtx)
  new.rnames = sapply(rnames, function(x) gsub('-', '_', x))
  rownames(mtx) <- new.rnames
  
  dt = reshape2::melt(as.matrix(mtx), value.name = 'count')
  rm(mtx)
  dt = dt[dt$count > 0, ]
  input_cds <- make_atac_cds(dt, binarize = T)
  rm(dt)
  input_cds <- detectGenes(input_cds)
  input_cds <- estimateSizeFactors(input_cds)
  
  if(reduction == 'tsne') {
    if(is.null(seurat.obj@reductions$tsne))
      seurat.obj <- RunTSNE(seurat.obj, dims = 1:npc)
    redu.coords = seurat.obj@reductions$tsne@cell.embeddings
  }
  if(reduction == 'umap') {
    if(is.null(seurat.obj@reductions$umap))
      seurat.object <- RunUMAP(seurat.object, dims = 1:npc)
    redu.coords = seurat.obj@reductions$umap@cell.embeddings
  }
  
  #make the cell id consistet
  
  cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = redu.coords)
  
  ## get connections
  
  conns <- run_cicero(cicero_cds, chr_sizes) 
  

  return(conns)
}



## query the resoltuion parameters given a seurat object and the number of clusters
## using binary seach
queryResolution4Seurat <- function(seurat.obj, k = 10, reduction = 'umap', npc = 20, 
                      min_resl = 0.1, max_resl = 1, max_iter = 15, doPCA = F){
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
    if(k1 == 10) stop('Please specify a much smaller min_res')
  }

  while(len2 < k){
    k2 = k2 + 1
    message('max_resl too small, trying to multiply it by 2')
    max_resl = max_resl * 2
    tmp.cluster2 <- FindClusters(seurat.obj, resolution = max_resl)@active.ident
    len2 = length(levels(tmp.cluster2))
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



## query the resoltuion parameters given a seurat object (transformed from cistopic object) and the number of clusters
## using binary seach
queryResolution4Topic <- function(seurat.obj, k = 10,  min_resl = 0.1, max_resl = 1, max_iter = 15){

  
  # skip find neighbors 
 
  tmp.cluster1 <- FindClusters(seurat.obj, resolution = min_resl, graph.name = 'snn')@active.ident
  tmp.cluster2 <- FindClusters(seurat.obj, resolution = max_resl, graph.name = 'snn')@active.ident
  
  
  

  len1 = length(levels(tmp.cluster1))
  len2 = length(levels(tmp.cluster2))

  k1 = k2 = 0
  while(len1 > k ){
   
    k1 = k1 + 1
    message('min_resl too large, trying to divided it by  2')
    min_resl = min_resl/2
    tmp.cluster1 <- FindClusters(seurat.obj, resolution = min_resl, graph.name = 'snn')@active.ident
    len1 = length(levels(tmp.cluster1))
    if(k1 == 10) stop('Please specify a much smaller min_res')
  }

  while(len2 < k){
    k2 = k2 + 1
    message('max_resl too small, trying to multiply it by 2')
    max_resl = max_resl * 2
    tmp.cluster2 <- FindClusters(seurat.obj, resolution = max_resl, graph.name = 'snn')@active.ident
    len2 = length(levels(tmp.cluster2))
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
    
    tmp.cluster <- FindClusters(seurat.obj, resolution = resl0, graph.name = 'snn')@active.ident
      
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
queryResolution4Topic = cmpfun(queryResolution4Topic)



# plot cicero gene activity score, given a gene set, cicero gene activity score and seurat atac object
plotCicero_ascore <- function(genes, cicero_ascore, seurat.atac, reduction = 'tsne'){
  ## genes should include gene_name and gene_id information 
  ## gene_id is the rownames of cicero_ascore, cell id is the colnames of cicero_ascore
  
  cicero_ascore = cicero_ascore[rownames(cicero_ascore) %in% genes$gene_id, ]
  
  if(nrow(cicero_ascore) == 0) stop('No activity score found for these genes!')
  genes = genes[gene_id %in% rownames(cicero_ascore), ]

  # add the gene activity score as a metadata feature
  if(any(colnames(seurat.atac) != colnames(cicero_ascore))) stop('Cell id not consistent!')

  for(gene_id0 in genes$gene_id){
    seurat.atac@meta.data[[gene_id0]] <- cicero_ascore[gene_id0, ]
  }

  tmp_plot = list()
  for(i in 1:nrow(genes)){
    tmp_plot[[i]] <- FeaturePlot(object = seurat.atac, reduction = reduction,  feature = genes$gene_id[i]) + 
    labs(title = genes$gene_name[i])
  } 
  
  pp <- ggarrange(plotlist = tmp_plot, ncol = 2, nrow = ceiling(length(tmp_plot)/2))
  return(pp)
}


basicCluster <- function(reduced.mtx, method = 'hclust', k = 5){
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



## do clustering using different #pcs and different methods, given an seurat.obj
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

     
     seurat.obj@meta.data[[clustLabelName]] = basicCluster(embd.mat, method, k)
  }


 
  return(seurat.obj)
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
  cisTopicObject <- runModels(cisTopicObject, topic = c(5, 10, 15, 20, 25, 30, 40, 50), seed = 987, nCores = nCores, 
    burnin = 120, iterations = 150, addModels = T)
  #cisTopicObject <- selectModel(cisTopicObject, keepBinarymatrix = F, keepModels = F)
  #cellassign <- t(modelMatSelection(cisTopicObject, 'cell', 'Probability'))
  return(cisTopicObject)
}


# the imput mtx is already filterd
run_scABC <- function(mtx, k = 5){
  weights = apply(mtx, 2, mean)
  landmarks = computeLandmarks(mtx, weights = weights, nCluster = k)
  labels = assign2landmarks(mtx, landmarks)
  return(labels)
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


##calculate rand index, given true labels and npc, and a seurat object
## return data.table recoords rand index for different conditions
calRand_gTrueLabelsAndNTOPIC <- function(cisTopic.obj, true.labels, ntopic = 20, 
                                      clust_methods = c('gcSNN', 'hclust', 'kmeans'), 
                                      reductions = c('topic', 'tsne', 'umap'), resolution = 0.1){
   

   rands = matrix(0, length(clust_methods), length(reductions))

   for(i in 1:nrow(rands)){

    if(clust_methods[i] == 'gcSNN'){
      for (j in 1:ncol(rands)){
      pred.labels <- gcSNN4CistopicObj(cisTopic.obj, ntopic = ntopic, clust_by = reductions[j], resolution = resolution, 
        k = length(unique(true.labels)))$cluster_label

      rands[i, j] = adjustedRandIndex(pred.labels, true.labels)
      }
    }else{
      for (j in 1:ncol(rands)){
        sele.cisTopic <- selectModel(cisTopic.obj, select = ntopic, 
                               keepBinaryMatrix = F, keepModels = F)
        sele.cisTopic <- runtSNE(sele.cisTopic, target='cell')
        sele.cisTopic <- runUmap(sele.cisTopic, target='cell')
        if(reductions[j] == 'topic') cell_topic <- t(modelMatSelection(sele.cisTopic, 'cell', 'Probability'))
        if(reductions[j] == 'tsne') cell_topic <- sele.cisTopic@dr$cell$tSNE
        if(reductions[j] == 'umap') cell_topic <- sele.cisTopic@dr$cell$Umap
        pred.labels <- basicCluster(cell_topic, method = clust_methods[i], k = length(unique(true.labels)))

        rands[i, j] = adjustedRandIndex(pred.labels, true.labels)
     }
   }
 }
   rownames(rands) = clust_methods
   colnames(rands) = reductions
   rands = reshape2::melt(rands, value.name = 'rand')
   rands$ntopic = ntopic

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


## compare clustering using different reductions, given npc and method  -- not used
compRand_reduction_gNPCAndMethod <- function(seurat.obj, clust_method = 'gcSNN', 
                                  reductions = c('pca', 'tsne', 'umap'),
                                  resolution = 0.2, k = 10, npc = 20, ...){

  if(is.null(seurat.obj@reductions$pca)){
    seurat.obj <- RunPCA(seurat.obj, npcs = npc, verbose = F)
  }

  for(reduction0 in reductions){
    seurat.obj <- clust4seurat(seurat.obj, npc = npc, method = clust_method, 
                                  reduction = reduction0, resolution = resolution, k = k,
                                  clustLabelName = paste(clust_method, '_', reduction0, '_npc', npc))

  }
  #embd.mtx = seurat.obj@reductions[['tsne']]@cell.embeddings

  # evaluat clustering, calculate corrected rand index
  rands = rep(0, length(reductions) - 1)
  for(i in 1:(length(rands))){
    #rands[i] = eval_cluster4mat(embd.mtx, 
    #seurat.obj@meta.data[[paste(clust_method, '_', reductions[i], '_npc', npc)]],
    #seurat.obj@meta.data[[paste(clust_method, '_', reductions[i+1], '_npc', npc)]])$rand
    rands[i] = adjustedRandIndex(seurat.obj@meta.data[[paste(clust_method, '_', reductions[i], '_npc', npc)]],
    seurat.obj@meta.data[[paste(clust_method, '_', reductions[i+1], '_npc', npc)]])
  }
  
  set.cols = brewer.pal(n = length(rands), name = 'Dark2')
  
  p <- barplot(rands, col = set.cols, ylab = 'Adjust Rand Index', ylim = c(0, 1),  ...)
  
  return(rands)
   
}

## compare clustering using different methods, given npc and reduction 
compRand_method_gNPCAndReduction <- function(seurat.obj, clust_methods = c('gcSNN', 'hclust', 'kmean'), 
                                  reduction = 'pca',
                                  resolution = 0.2, k = 10, npc = 20, plotDR = FALSE, ...){
  if(is.null(seurat.obj@reductions$pca)){
    seurat.obj <- RunPCA(seurat.obj, npcs = npc, verbose = F)
  } 
  i = 0
  pp_tsne = pp_umap = list()
  for(method0 in clust_methods){
    i = i + 1
    seurat.obj <- clust4seurat(seurat.obj, npc = npc, method = method0, 
                                  reduction = reduction, resolution = resolution, k = k,
                                  clustLabelName = paste(method0, '_', reduction, '_npc', npc))
    if(plotDR) {
       seurat.obj = RunTSNE(seurat.obj, dims = 1:npc)
      pp_tsne[[i]] = DimPlot(seurat.obj, reduction = 'tsne', group.by = paste(method0, '_', reduction, '_npc', npc))
       seurat.obj = RunUMAP(seurat.obj, dims = 1:npc, verbose = F)
      pp_umap[[i]] = DimPlot(seurat.obj, reduction = 'umap', group.by = paste(method0, '_', reduction, '_npc', npc))

    }
  }
  #embd.mtx = seurat.obj@reductions[[reduction]]@cell.embeddings

  # evaluat clustering, calculate corrected rand index
  rands = rep(0, length(clust_methods) - 1)
  for(i in 1:(length(rands))){
    #rands[i] = eval_cluster4mat(embd.mtx, 
    #seurat.obj@meta.data[[paste(clust_methods[i], '_', reduction, '_npc', npc)]],
    #seurat.obj@meta.data[[paste(clust_methods[i + 1], '_', reduction, '_npc', npc)]])$rand
    rands[i] = adjustedRandIndex(seurat.obj@meta.data[[paste(clust_methods[i], '_', reduction, '_npc', npc)]],
    seurat.obj@meta.data[[paste(clust_methods[i + 1], '_', reduction, '_npc', npc)]])

  }
  
  res = list('rand' = rands)
  if(plotDR){
    res$plots_tsne = pp_tsne
    res$plots_umap = pp_umap
  }
  return(res)

}




## compare clustering rand index using different npc, given a method and reduction
compRand_ntopic_gReductionAndMethod <- function(cisTopic.obj, clust_method = 'gcSNN', reduction = 'topic',
                                resolution = 0.2, k = 10, ntopics = c(20, 30, 40, 50), ...){

  
  snn.res = list()
  i=1
  for(ntopic0 in ntopics){
    if(clust_method == 'gcSNN') {
      snn.res[[i]] <- gcSNN4CistopicObj(cisTopic.obj, ntopic0, reduction, resolution = resolution, k = k)
      }else{
        
        

        sele.cisTopic <- selectModel(cisTopic.obj, select = ntopic0, 
                               keepBinaryMatrix = F, keepModels = F)
        sele.cisTopic <- runtSNE(sele.cisTopic, target='cell')
        sele.cisTopic <- runUmap(sele.cisTopic, target='cell')
        if(reduction == 'topic') cell_topic <- t(modelMatSelection(sele.cisTopic, 'cell', 'Probability'))
        if(reduction == 'tsne') cell_topic <- sele.cisTopic@dr$cell$tSNE
        if(reduction == 'umap') cell_topic <- sele.cisTopic@dr$cell$Umap
        snn.res[[i]] <- basicCluster(cell_topic, method = clust_method, k = k)

      }
    
    i = i+1
    
  }

  # evaluat clustering, calculate corrected rand index
  rands = rep(0, length(ntopics) - 1)
  for(i in 1:length(rands)){
    if(clust_method == 'gcSNN') {
        rands[i] = adjustedRandIndex(snn.res[[i]]$cluster_label,
                                snn.res[[i+1]]$cluster_label)
    }else{
        rands[i] = adjustedRandIndex(snn.res[[i]],
                                snn.res[[i+1]])
    }
  }
  
  set.cols = brewer.pal(n = length(rands), name = 'Dark2')
  
  p <- barplot(rands, col = set.cols, ylab = 'Adjust Rand Index', ylim = c(0, 1), ...)
  
  return(rands)

}


## compare clustering rand index using different npc, given a method and reduction
compRand_method_gNTOPICAndReduction <- function(cisTopic.obj, clust_methods = c('gcSNN', 'hclust', 'kmeans'), reduction = 'topic',
                                resolution = 0.2, k = 10, ntopic = 30, ...){

  
  snn.res = list()
  i=1
  for(method0 in clust_methods){
    if(clust_methods[i] == 'gcSNN') {
      snn.res[[i]] <- gcSNN4CistopicObj(cisTopic.obj, ntopic, reduction, resolution = resolution, k = k)
      }else{
        
        

        sele.cisTopic <- selectModel(cisTopic.obj, select = ntopic, 
                               keepBinaryMatrix = F, keepModels = F)
        sele.cisTopic <- runtSNE(sele.cisTopic, target='cell')
        sele.cisTopic <- runUmap(sele.cisTopic, target='cell')
        if(reduction == 'topic') cell_topic <- t(modelMatSelection(sele.cisTopic, 'cell', 'Probability'))
        if(reduction == 'tsne') cell_topic <- sele.cisTopic@dr$cell$tSNE
        if(reduction == 'umap') cell_topic <- sele.cisTopic@dr$cell$Umap
        snn.res[[i]] <- basicCluster(cell_topic, method = method0, k = k)

      }
    
    i = i+1
    
  }


  # evaluat clustering, calculate corrected rand index
  rands = rep(0, length(clust_methods) - 1)
  for(i in 1:length(rands)){
    if(clust_methods[i] == 'gcSNN') {
        rands[i] = adjustedRandIndex(snn.res[[i]]$cluster_label,
                                snn.res[[i+1]]$cluster_label)
    }else{
        rands[i] = adjustedRandIndex(snn.res[[i]],
                                snn.res[[i+1]])
    }
  }
  
  set.cols = brewer.pal(n = length(rands), name = 'Dark2')
  
  p <- barplot(rands, col = set.cols, ylab = 'Adjust Rand Index', ylim = c(0, 1), ...)
  
  return(rands)

}


## plot cisTopic cluster on dimension reduction plot
## transfer to seurat object to have a similar color scheme
plotDR4cisTopic_bySeurat <- function(cisTopic.obj, ntopic = 20, cistopic.cl.label){
  sele.cisTopic <- selectModel(cisTopic.obj, select = ntopic, 
                               keepBinaryMatrix = F, keepModels = F)
  sele.cisTopic <- runtSNE(sele.cisTopic, target='cell')
  sele.cisTopic <- runUmap(sele.cisTopic, target='cell')
  colnames(sele.cisTopic@dr$cell$tSNE) = c('tSNE_1', 'tSNE_2')
  colnames(sele.cisTopic@dr$cell$Umap) = c('UMAP_1', 'UMAP_2')
  


  embd.mtx = modelMatSelection(sele.cisTopic, 'cell', 'Probability')

  seurat.obj <- CreateSeuratObject(embd.mtx)
  seurat.obj@meta.data[['cistopic.cl.label']] = cistopic.cl.label
  seurat.obj = ScaleData(seurat.obj, do.scale = F, do.center = F, vars.to.regress = NULL)
  seurat.obj = RunPCA(seurat.obj, npc = floor(ntopic/2), features = rownames(seurat.obj), verbose = F)


  seurat.obj = RunTSNE(seurat.obj, dims = 1:floor(ntopic/2))
  #replace the tsne coordinate
  seurat.obj@reductions$tsne@cell.embeddings = sele.cisTopic@dr$cell$tSNE
  pp_tsne = DimPlot(seurat.obj, reduction = 'tsne', group.by = 'cistopic.cl.label')

  
  seurat.obj = RunUMAP(seurat.obj, dims = 1:floor(ntopic/2), verbose = F)
  seurat.obj@reductions$umap@cell.embeddings = sele.cisTopic@dr$cell$Umap
  #replace the umap coordinate
  pp_umap = DimPlot(seurat.obj, reduction = 'umap', group.by = 'cistopic.cl.label')



  return(list('plot_tsne' = pp_tsne, 'plot_umap' = pp_umap, 'seurat.obj' = seurat.obj))
}



