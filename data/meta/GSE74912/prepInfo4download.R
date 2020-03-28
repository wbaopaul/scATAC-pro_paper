library(data.table)

dd = fread('SraRunInfo.csv')
ss = fread('sra_result.csv')

dd = subset(dd, select = c('Run', 'Sample', 'SampleName'))
ss = subset(ss, select = c('Experiment Accession', 'Experiment Title', 
                           'Sample Accession'))
names(ss) = c('exp_acc', 'title', 'Sample')
setkey(ss, Sample)
dd[, 'title' := ss[J(dd$Sample)]$title]
dd[, title := unlist(strsplit(title, ';'))[1], by = title]
dd[, title := unlist(strsplit(title, ':'))[2], by = title]
dd[, 'exp_acc' := ss[J(dd$Sample)]$exp_acc]

dd[, 'donor' := unlist(strsplit(title, '-'))[1], by = title]
dd[, 'ctype' := unlist(strsplit(title, '-'))[2], by = title]
dd[, donor := gsub('Donor', '', donor), by = donor]
dd[, ctype := gsub('1A|1B', 'HSC', ctype), by = ctype]
dd[, ctype := gsub('2A|2B', 'MPP', ctype), by = ctype]
dd[, ctype := gsub('3A|3B', 'LMPP', ctype), by = ctype]
dd[, ctype := gsub('4A|4B', 'CMP', ctype), by = ctype]
dd[, ctype := gsub('5A|5B', 'GMP', ctype), by = ctype]
dd[, ctype := gsub('6A|6B', 'MEP', ctype), by = ctype]
dd[, ctype := gsub('7A|7B', 'Mono', ctype), by = ctype]
dd[, ctype := gsub('Nkcell', 'NK', ctype), by = ctype]

dd = dd[ctype != 'CD34']

dd = dd[!grepl(dd$title, pattern = 'SU'),]  ## remove cancer cells
write.table(dd, file = 'srrInfo4Download.txt', sep = '\t', row.names = F, quote = F)
