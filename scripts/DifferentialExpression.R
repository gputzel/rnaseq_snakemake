library(tidyverse)
library(DESeq2)
library(ggp.rnaseq.human)

comparison <- snakemake@wildcards[['comparison']]

params <- snakemake@config[['comparisons']][[comparison]]
bad.samples <- params[['bad-samples']]
new.design <- as.formula(params[['design']])
prefilter.string <- params[['prefilter']]
sample.subset <- params[['subset']]
de.attribute <- params[['attribute']]
de.group.1 <- params[['group1_plus']]
de.group.2 <- params[['group2_minus']]

dds <- readRDS(snakemake@input[['dds']])

coldata <- colData(dds)
coldata$sizeFactor <- NULL #Get rid of size factors estimated using all samples
countdata <- counts(dds,normalize=F)

dds.new <- DESeqDataSetFromMatrix(countData = countdata,colData=coldata,design=new.design)

dds.no.bad.samples <- dds.new[,!(colnames(dds.new) %in% names(bad.samples))]

sample.selection <- with(colData(dds.no.bad.samples),eval(parse(text=sample.subset)))

dds.subset <- dds.no.bad.samples[,sample.selection]

for(variable in all.vars(design(dds.subset))){
    cat(variable,'\n')
    colData(dds.subset)[,variable] <- droplevels(colData(dds.subset)[,variable])
}

eval(parse(text=paste0('prefilter <- ',prefilter.string)))

dds.filtered <- prefilter(dds.subset)

dds.filtered <- DESeq(dds.filtered)

saveRDS(dds.filtered,snakemake@output[['dds']])

res <- ggp.rnaseq.human::ggp.results(dds.filtered,c(de.attribute,de.group.1,de.group.2))

saveRDS(res,snakemake@output[['rds']])
