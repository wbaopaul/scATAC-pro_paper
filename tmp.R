rand.clean = readRDS(file = 'data/intermediate/erythropoiesis_clean/eval_binary_rand_erythropoiesis_clean_PC30.rds')
rand.noisy2 = readRDS(file = 'data/intermediate/erythropoiesis_noisy_p2/eval_binary_rand_erythropoiesis_noisy_p2_PC30.rds')
rand.noisy4 = readRDS(file = 'data/intermediate/erythropoiesis_noisy_p4/eval_binary_rand_erythropoiesis_noisy_p4_PC30.rds')

## louvain ####
nvegs = c(5000, 10000, 15000, 20000, 25000)
postscript('output/figures/erythropoiesis_noisy_p4/randIndex_compare_bi_nonbi_louvain_erythropoiesis_noisy_p4.eps',
           height = 6, width = 6)
plot(nvegs, rand.noisy4$rand.lv[, 1], pch = 17, ylim = c(0, 1),
     ylab = 'Adjusted Rand', xlab = '#variable features',
     col = 2)
points(nvegs, rand.noisy4$rand.lv[, 2], pch = 19, col = 1)
legend(x = 19000, y = 0.4, legend = c('non-binary', 'binary'),
       col = c(2, 1), pch = c(17, 19), bty = 'o')
dev.off()

postscript('output/figures/erythropoiesis_noisy_p2/randIndex_compare_bi_nonbi_louvain_erythropoiesis_noisy_p2.eps',
           height = 6, width = 6)
plot(nvegs, rand.noisy2$rand.lv[, 1], pch = 17, ylim = c(0, 1),
     ylab = 'Adjusted Rand', xlab = '#variable features',
     col = 2)
points(nvegs, rand.noisy2$rand.lv[, 2], pch = 19, col = 1)
legend(x = 19000, y = 0.4, legend = c('non-binary', 'binary'),
       col = c(2, 1), pch = c(17, 19), bty = 'o')
dev.off()

postscript('output/figures/erythropoiesis_clean/randIndex_compare_bi_nonbi_louvain_erythropoiesis_clean.eps',
           height = 6, width = 6)
plot(nvegs, rand.clean$rand.lv[, 1], pch = 17, ylim = c(0, 1),
     ylab = 'Adjusted Rand', xlab = '#variable features',
     col = 2)
points(nvegs, rand.clean$rand.lv[, 2], pch = 19, col = 1)
legend(x = 19000, y = 0.4, legend = c('non-binary', 'binary'),
       col = c(2, 1), pch = c(17, 19), bty = 'o')
dev.off()

## kmeans ####

postscript('output/figures/erythropoiesis_noisy_p4/randIndex_compare_bi_nonbi_kmeans_erythropoiesis_noisy_p4.eps',
           height = 6, width = 6)
plot(nvegs, rand.noisy4$rand.kmeans[, 1], pch = 17, ylim = c(0, 1),
     ylab = 'Adjusted Rand', xlab = '#variable features',
     col = 2)
points(nvegs, rand.noisy4$rand.kmeans[, 2], pch = 19, col = 1)
legend(x = 19000, y = 0.4, legend = c('non-binary', 'binary'),
       col = c(2, 1), pch = c(17, 19), bty = 'o')
dev.off()

postscript('output/figures/erythropoiesis_noisy_p2/randIndex_compare_bi_nonbi_kmeans_erythropoiesis_noisy_p2.eps',
           height = 6, width = 6)
plot(nvegs, rand.noisy2$rand.kmeans[, 1], pch = 17, ylim = c(0, 1),
     ylab = 'Adjusted Rand', xlab = '#variable features',
     col = 2)
points(nvegs, rand.noisy2$rand.kmeans[, 2], pch = 19, col = 1)
legend(x = 19000, y = 0.4, legend = c('non-binary', 'binary'),
       col = c(2, 1), pch = c(17, 19), bty = 'o')
dev.off()

postscript('output/figures/erythropoiesis_clean/randIndex_compare_bi_nonbi_kmeans_erythropoiesis_clean.eps',
           height = 6, width = 6)
plot(nvegs, rand.clean$rand.kmeans[, 1], pch = 17, ylim = c(0, 1),
     ylab = 'Adjusted Rand', xlab = '#variable features',
     col = 2)
points(nvegs, rand.clean$rand.kmeans[, 2], pch = 19, col = 1)
legend(x = 19000, y = 0.4, legend = c('non-binary', 'binary'),
       col = c(2, 1), pch = c(17, 19), bty = 'o')
dev.off()

