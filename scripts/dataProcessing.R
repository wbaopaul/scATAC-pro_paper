library(data.table)
library(magrittr)
library(readr)
library(tidyr)
library(Seurat)

## read raw data

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


readRawData <- function(data_name, minFrac_in_cell = 0.01, min_depth = 1000){
  if(data_name == 'HSC_GSE96769'){
    # load scATAC count data
    dt = fread("data/raw/HSC_GSE96769/GSE96769_scATACseq_counts.txt", skip = 1)
    cname = readLines('data/raw/HSC_GSE96769/GSE96769_scATACseq_counts.txt', n = 1)
    cname = unlist(strsplit(cname, ';'))
    
    names(dt) = c('peakID', 'cellID', 'count')
    
    # load peak annotation
    peaks = fread('data/raw/HSC_GSE96769/GSE96769_PeakFile_20160207.bed',
                  select = 1:3)
    names(peaks) = c('chr', 'start', 'end')
    
  
    dt[, 'cell.freq' := as.double(sum(count)), by = cellID]
    cellID.filterd = unique(dt[cell.freq < min_depth]$cellID)
    if(!is.null(cellID.filterd)) dt = dt[!cellID %in% cellID.filterd]
    
    dt[, 'peak.freq' := .N, by = peakID]
    peakID.filterd = unique(dt[peak.freq < minFrac_in_cell * length(cname)]$peakID)
    if(!is.null(peakID.filterd)) dt = dt[!peakID %in% peakID.filterd]
    
    
    dt[, c('cell.freq', 'peak.freq') := NULL]
    atac.mtx = dcast(dt,  peakID ~ cellID, value.var = "count")
    rm(dt)
    rownames(atac.mtx) = atac.mtx$peakID
    atac.mtx[, 'peakID' := NULL]
    atac.mtx[is.na(atac.mtx)] = 0L
    colnames(atac.mtx) = cname[as.integer(colnames(atac.mtx))]
    
    
    pos.peaks = peaks[as.integer(rownames(atac.mtx))]
    rownames(atac.mtx) = pos.peaks[, 'pos' := paste(chr, start, end, sep = '-')]$pos
    
    feature.name = rownames(atac.mtx)
    atac.mtx = as.matrix(atac.mtx)
    rownames(atac.mtx) = feature.name
    atac.mtx = as(atac.mtx, 'sparseMatrix')
  }
  
  if(data_name == 'PMBC10k'){
    dirt = './data/raw/PBMC10k/filtered_peak_bc_matrix/'
    atac.mtx = read10X_ATAC(dirt)
    depth.cell = Matrix::colSums(atac.mtx)
    atac.mtx = atac.mtx[, depth.cell > min_depth]
    frac.in.cell = Matrix::rowSums(atac.mtx > 0)
    atac.mtx = atac.mtx[frac.in.cell > minFrac_in_cell, ]
  }
  return(atac.mtx)
}


filterMat <- function(atac.mtx, minFrac_in_cell = 0.01, min_depth = 1000){
  depth.cell = Matrix::colSums(atac.mtx)
  atac.mtx = atac.mtx[, depth.cell > min_depth]
  frac.in.cell = Matrix::rowSums(atac.mtx > 0)
  atac.mtx = atac.mtx[frac.in.cell > minFrac_in_cell, ]
  return(atac.mtx)
}

## tf-idf normalization
atac_tfidf = function(atac_matrix, site_frequency_threshold=0.03) {
  num_cells_ncounted = Matrix::rowSums(atac_matrix)
  threshold = ncol(atac_matrix) * site_frequency_threshold
  
  ncounts = atac_matrix[num_cells_ncounted >= threshold,]
  
  ## Normalize the data with TF-IDF
  nfreqs = t(t(ncounts) / Matrix::colSums(ncounts))
  tf_idf_counts = nfreqs * log(1 + ncol(ncounts) / Matrix::rowSums(ncounts))
  
  return(list('tfidf_matrix' = tf_idf_counts, 'filtered_matrix' = ncounts))
}
