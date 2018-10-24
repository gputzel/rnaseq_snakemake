library(ggp.rnaseq)
library(DESeq2)

dds <- readRDS(snakemake@input[['rds']])

dds <- DESeq(dds)

ggp.rnaseq::write.normalized.counts(dds,filename=snakemake@output[['counts']],sample.attributes=c('genotype','cell.type','mouse'),overwrite=T)
