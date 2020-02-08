library(splatter)
library(scater)

generateSCAtac_splatter <- function(integer.mtx, group.prob = c(0.3, 0.3, 0.4),
                                    de.prob = c(0.3, 0.3, 0.3), params = NULL,
                                    ncells = 5000, npeaks = 20000, isNorm = T){
  if(is.null(params)){
    params = newSplatParams()
    params <- setParams(params, update = list(lib.norm = isNorm))
    params = splatEstimate(integer.mtx, params = params)
  }
  params <- setParams(params, update = list(nGenes = npeaks, batchCells = ncells,
                                            lib.norm = isNorm))
  simu = splatSimulate(group.prob = group.prob, method = 'groups', params = params,
                       de.prob = de.prob)
  #simu = normalize(simu)
  #plotTSNE(simu, colour_by = 'Group')
  return(list('mtx' = simu@assays$data$counts, 
              'cluster_labels' = simu@colData[, c('Cell', 'Group')],
              'params' = params))
}


## generate synthetic data with noise (download from with some minor modification)
#-----------
# Parameters
#-----------
# n_cells Number of cells from each group to be simulated; either a number (all the same) or a vector
# of length(which_celltypes)

# which_celltypes Number of groups to partition cells into (even groups); must evenly divide into n_cells
# n_frags_per_cell number of fragments in peaks to be simulated per single cell

# rate_noise number between 0 (perfect downsample) and 1 (nonsense) for noise

# seed random parameter for setting the seed
# shuffle Randomly order the resulting cells and peaks
simulate_scatac <- function(bulk, n_cells, which_celltypes, n_frags_per_cell = 1000, 
                            rate_noise = 0, seed = 100, shuffle = FALSE){
  
  # Reproducibility
  set.seed(seed)
  which_celltypes <- sort(which_celltypes)
  stopifnot(rate_noise < 1) 
  stopifnot(n_frags_per_cell > 100)
  n_peaks <- dim(bulk)[1]
  #--
  # Set up cell labels
  #--
  
  if(length(n_cells) > 1){
    stopifnot(length(which_celltypes) == length(n_cells))
    
    # Generate cell labels
    cell_labels <- sapply(1:length(which_celltypes), function(i){
      rep(which_celltypes[i], n_cells[i])
    }) %>% unlist() %>% sort()
    
  } else {
    #     n_groups <- length(which_celltypes)
    #     cell_labels <- sort(rep(which_celltypes, n_cells*n_groups))
    cell_labels <- sort(rep(which_celltypes, n_cells))
  }
  final_names <- paste0(cell_labels, "_", as.character(1:length(cell_labels)))
  
  
  #-------------------
  # Simulate true data
  #-------------------
  
  # Generate cell-type specific peaks
  lapply(which_celltypes, function(celltype){
    
    # Apply different rates per cell depending on group label for generating cell-type specific peaks
    n_cells_this_celltype <- sum(cell_labels == celltype)
    counts_celltype <- bulk[,celltype]
    
    # Define probabilities
    #                        Prob observting frag                Total number of fragments epxpected; the 0.5s are for two alleles that will be simulated/added later
    prob_per_peaks <- counts_celltype/sum(counts_celltype) * (n_frags_per_cell*0.5 * (1-rate_noise)) + ((rate_noise*n_frags_per_cell)/n_peaks*0.5) 
    
    # Cap probabilities at something sensible
    prob_per_peaks <- ifelse(prob_per_peaks > 0.9, 0.9, prob_per_peaks)
    
    # Represent the two haplotypes as two random draws
    mat1 <- (matrix(rbinom(n_peaks*n_cells_this_celltype, size = 1, prob = prob_per_peaks),
                    ncol = n_cells_this_celltype, byrow = FALSE) )
    mat2 <- (matrix(rbinom(n_peaks*n_cells_this_celltype, size = 1, prob = prob_per_peaks),
                    ncol = n_cells_this_celltype, byrow = FALSE) )
    
    mat <- mat1 + mat2
    Matrix(mat)
  }) %>% do.call(what = "cbind") -> sparse_matrix
  
  colnames(sparse_matrix) <- final_names
  #peaknames = paste(peaks$V1,peaks$V2,peaks$V3,sep = "_")
  #rownames(sparse_matrix) <- peaknames
  rownames(sparse_matrix) <- rownames(bulk.mtx)
  sparse_matrix
}
