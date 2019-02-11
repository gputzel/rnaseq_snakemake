library(ggp.rnaseq.human)
library(DESeq2)

dds <- readRDS(snakemake@input[['rds']])

dds <- DESeq(dds)

ggp.rnaseq.human::write.normalized.counts(dds,filename=snakemake@output[['counts']],sample.attributes=c('Group','Patient','Site'),overwrite=T)
