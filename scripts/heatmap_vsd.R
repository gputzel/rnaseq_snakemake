library(DESeq2)

heatmap <- snakemake@wildcards[['heatmap']]
params <- snakemake@config[['heatmaps']][[heatmap]]
prefilter.string <- params[['prefilter']]
sample.subset <- params[['subset']]
bad.samples <- params[['bad-samples']]
new.design <- as.formula(params[['design']])

dds <- readRDS(snakemake@input[['dds']])

coldata <- colData(dds)
coldata$sizeFactor <- NULL #Get rid of size factors estimated using all samples
countdata <- counts(dds,normalize=F)

dds.new <- DESeqDataSetFromMatrix(countData = countdata,colData=coldata,design=new.design)

dds.no.bad.samples <- dds.new[,!(colnames(dds.new) %in% names(bad.samples))]
sample.selection <- with(colData(dds.no.bad.samples),eval(parse(text=sample.subset)))
dds.subset <- dds.no.bad.samples[,sample.selection]

eval(parse(text=paste0('prefilter <- ',prefilter.string)))
dds.filtered <- prefilter(dds.subset)

dds.filtered <- DESeq(dds.filtered)
vsd <- varianceStabilizingTransformation(dds.filtered)
saveRDS(vsd,file=snakemake@output[['vsd']])