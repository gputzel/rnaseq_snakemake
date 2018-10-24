library(ggp.rnaseq)
library(ggp.genesets)
library(tidyverse)
library(fgsea)
library(DESeq2)

fgsea.params <- snakemake@config[['fgsea']]

res <- readRDS(snakemake@input[['rds']])

res %>%
    dplyr::select(ENTREZID,stat) %>%
    na.omit() %>%
    deframe -> ranks

all.go.names %>% enframe(value="Pathway",name="GOID") -> go.names

all.go.unique <- lapply(all.go,unique)

res.fgsea <- fgsea(pathways=all.go.unique,stats=ranks,
                   minSize=fgsea.params[['minSize']],
                   maxSize=fgsea.params[['maxSize']],
                   nperm=fgsea.params[['nperm']])

inner_join(go.names,res.fgsea,by=c("GOID"="pathway")) -> res.fgsea.with.names

saveRDS(res.fgsea.with.names,snakemake@output[['fgsea']])
