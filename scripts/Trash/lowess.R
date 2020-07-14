library(bbmle)
library(emdbook)
library(limma)

sysargv = commandArgs(trailingOnly=TRUE)
BAD = sysargv[1]
sp = sysargv[2]
weights_file = paste0('~/PARAMETERS/weights_BAD=', BAD, '.tsv')

weights_df = read.table(weights_file, header=TRUE)

weights_df$lowess_r_weights = weightedLowess(as.numeric(row.names(weights_df)), weights_df$alpha, weights=weights_df$snps, span=as.numeric(sp))$fitted

write.table(weights_df, file=paste0('~/PARAMETERS/r_weights_BAD=', BAD, '.tsv'), sep='\t', row.names=FALSE)

