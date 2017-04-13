#!/usr/bin/env Rscript
# exctract sample IDs from Illumina FinalReport file in Genomelist format

library(data.table)
args = commandArgs(trailingOnly=TRUE)
write(paste(Sys.time(),' exctracting IDs from ',args[1]),stderr())
ids <- fread(args[1],select=2,col.names='Id',showProgress=F)[,1,with=F]
write.table(data.table(File=args[1],unique(ids)),stdout(),row.names=F,col.names=F,quote=F,sep='\t')


