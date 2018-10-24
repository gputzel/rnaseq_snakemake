library(ggp.rnaseq) #Commit 7d9b299fdfd10414a428bcb5d84e5b6d2b3b3419
library(tidyverse)
library(DESeq2)

counts <- read.table(snakemake@input[['counts']], sep = "\t", 
            header = TRUE, row.names = 1)

colnames(counts) %>%
    enframe(name='X',value = 'filename') %>%
    mutate(filename=gsub("X[.]Volumes[.]IBD[.]Putzel[.]Pubshare_backup[.]Artis_Lab[.]AnneLaure_Flamar[.]RNAseq[.]Tph1_KO_MLN[.]BAM[.]","",filename)) %>%
    mutate(filename=gsub("_PF[.]aligned[.]sorted[.]bam","",filename)) %>% 
    mutate(filename=gsub("^13i","Tph1_cKO_iILC2_13",filename)) %>% #Fix names that got screwed up
    mutate(filename=gsub("^13n","Tph1_cKO_mILC2_13",filename)) %>%
    mutate(filename=gsub("Tph1_cKO","KO",filename)) %>%
    mutate(filename=gsub("IL7R","cre",filename)) %>%
    extract(col = filename,into=c('genotype','cell.type','mouse'),regex="^([^_]*)_([^_]*)_([^_]*)",remove = F) %>%
    mutate(genotype = factor(genotype,levels=c('cre','KO'))) %>%
    mutate(cell.type = gsub("mILC2","nILC2",cell.type)) %>%
    mutate(cell.type = factor(cell.type,levels=c('nILC2','iILC2'))) %>%
    column_to_rownames(var = 'filename') %>% select(-X) -> coldata

coldata$mapping.rate <- import.mapping.rates()$mapping.rate

colnames(counts) <- rownames(coldata)

dds <- DESeqDataSetFromMatrix(countData=counts,
                              colData=coldata,
                              design=~genotype + cell.type)
dds <- DESeq(dds)

saveRDS(dds,file=snakemake@output[['rds']])
