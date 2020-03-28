library(data.table)
library(magrittr)
library(readr)
library(tidyr)
library(Seurat)
library(R.utils)
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


getRealData <- function(data_name, minFrac_in_cell = 0, min_depth = 1000,
                         max_depth = 50000){
  if(data_name == 'HSC_GSE96769'){
    # download scATAC count data
    count_file = 'data/raw/HSC_GSE96769/GSE96769_scATACseq_counts.txt.gz'
    peak_file = 'data/raw/HSC_GSE96769/GSE96769_PeakFile_20160207.bed.gz'
    download.file('https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE96769&format=file&file=GSE96769%5FscATACseq%5Fcounts%2Etxt%2Egz', 
                  count_file)
    download.file('https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE96769&format=file&file=GSE96769%5FPeakFile%5F20160207%2Ebed%2Egz', 
                  peak_file)
    
    dt = fread(count_file, skip = 1)
    cname = readLines(count_file, n = 1)
    cname = unlist(strsplit(cname, ';'))
    
    names(dt) = c('peakID', 'cellID', 'count')
    
    # download peak annotation
    peaks = fread(peak_file, select = 1:3)
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
  
  if(data_name == 'PBMC10k'){
    ## cellranger-atac processed matrix
    dirt = './data/raw/PBMC10k/filtered_peak_bc_matrix/'
    pbmc_url = 'http://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5'
    pbmc_file = 'data/raw/PBMC10k/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5'
    dir.create('data/raw/PBMC10k/', showWarnings = F)
    download.file(pbmc_url, pbmc_file)
    atac.mtx = Read10X_h5(pbmc_file)
    atac.mtx = filterMat(atac.mtx, min_depth = min_depth, 
                         max_depth = max_depth, minFrac_in_cell = minFrac_in_cell)
  }
  return(atac.mtx)
}

read_mtx_scATACpro <- function(mtx_path){
  #mtx_path <- paste0(dirt, "matrix.mtx")
  mtx.dir = dirname(mtx_path)
  feature_path <- paste0(mtx.dir, "/features.txt")
  barcode_path <- paste0(mtx.dir, "/barcodes.txt")
  
  
  features <- fread(feature_path, header = F)
  barcodes <- fread(barcode_path, header = F)
  
  mtx <-  Matrix::readMM(mtx_path) %>%
    magrittr::set_rownames(features$V1)%>%
    magrittr::set_colnames(barcodes$V1)
  
  return(mtx)
}


filterMat <- function(atac.mtx, minFrac_in_cell = 0.01, min_depth = 1000,
                      max_depth = 100000){
  depth.cell = Matrix::colSums(atac.mtx)
  atac.mtx = atac.mtx[, depth.cell > min_depth & depth.cell < max_depth]
  frac.in.cell = Matrix::rowMeans(atac.mtx > 0)
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
