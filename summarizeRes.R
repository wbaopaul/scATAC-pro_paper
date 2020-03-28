source('scripts/evalClusterMethods.R')
library(RColorBrewer)

## single run ####
dtype = 'HSC_GSE96769'  ## resample10k_hsc/HSC_GSE96769/GSE74912_clean/GSE74912_noisy_q2/GSE74912_noisy_q2_q4
mtx = readRDS(paste0('data/filtered_mat/', dtype, '_filtered_mat.rds'))
cnames = colnames(mtx)
true.label = sapply(cnames, function(x) unlist(strsplit(x, '_'))[1])
if(dtype == 'HSC_GSE96769') true.label = readRDS('output/intermediate/HSC_GSE96769/cell_type.rds')

methods = c('seurat', 'seurat_correct',   'cisTopic', 'scABC', 'chromVAR', 'LSI', 'SCRAT')
res = sumMethods4default(true.label, dtype, methods)
res = updateSummary(true.label, dtype, 'chromVAR', 
                    res$rands, res$pred.labels)
rands = res$rands
pred.labels = res$pred.labels

postscript(file = paste0('output/figures/', dtype, '/rand4', dtype, '.eps'),
           width = 8, height = 6)

set.cols = brewer.pal(n = length(rands), name = 'Dark2')

barplot(do.call('c', rands), col = set.cols, ylab = 'Adjust Rand', ylim = c(0, 1) )

dev.off()

## comapre time consuming ####
tt = fread('output/intermediate/resample10k_hsc/time4resample10k_hsc.txt')
methods = c('seurat', 'seurat_correct', 'cisTopic', 'scABC', 'chromVAR', 'LSI', 'SCRAT')
setkey(tt, V1)
tt = tt[methods, ]

postscript(file = paste0('output/figures/resample10k_hsc/time4resample10k_hsc.eps'),
           width = 8, height = 6)
set.cols = brewer.pal(n = length(tt$V2), name = 'Dark2')
barplot(tt$V2, col = set.cols, ylab = 'Computational Time (s)', names.arg = tt$V1 )
dev.off()

## run with different composition -- Fig S2C ####
dtype = 'resample10k_hsc'
rands = list()
methods = c('seurat', 'seurat_correct', 'cisTopic', 'scABC', 'chromVAR', 'LSI', 'SCRAT')

for(method0 in methods){
  fname = paste0('output/intermediate/resample10k_hsc/randIndexs/rand4', 
                 method0, '_resample10k_hsc.txt')
  if(!file.exists(fname)) next
  rands[[method0]] = fread(fname)
}


sRands = do.call('cbind', rands)
set.cols = brewer.pal(n = length(rands), name = 'Dark2')
postscript(file = paste0('output/figures/', dtype, '/rand4', dtype, '_diffComposition.eps'),
           width = 8, height = 6)
boxplot(sRands, names = names(rands), col = set.cols, ylab = 'Adjust Rand')
dev.off()



## check the composition
set.seed(2019)
props = rdirichlet(5, rep(3, 13))
groups = unique(true.label)

par(mar=c(6, 4.1, 4.1, 2.1))
boxplot(props, ylab = 'proportion', xlab = '', xaxt = 'n')
axis(1, labels = FALSE)
text(x = 1:13, y = par("usr")[3] - 0.01, srt = 45, adj = 1,
     labels = groups, xpd = TRUE)

set.cols = brewer.pal(n = ncol(props), name = 'Dark2')


groups[groups == 'Erythroblast'] = 'Ery'
postscript(file = paste0('output/figures/composition_', dtype, '_diffComposition.eps'),
           width = 10, height = 8)
par(xpd = T, mar = par()$mar + c(4,0,0,0))
barplot(t(props), names.arg = paste0('run', 1:5), col = set.cols,
        legend.text = groups, args.legend = 
          list(cex = 0.7, bty = 'n', horiz = T, text.width = 0.25, 
               x = 'bottom', x.intersp = 0.5, inset = c(0, -0.3)),
        ylab = 'cell type composition')
dev.off()
