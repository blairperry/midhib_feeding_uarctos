

library(tidyverse)
library(DESeq2)
library(IHW)
library(patchwork)

# Read in and reformat raw count table ------------------------------------

raw <- read_tsv('data/3_featureCounts/postDexExperiment_cellCultureData_02.10.23.txt',comment = '#')

# Clean sample names
names(raw) <- str_remove_all(names(raw),'/data/kelley/projects/bear_dextrose_reversal_newAnalysis_bwp/2_star/star_mapped/|Aligned.sortedByCoord.out.bam')

raw.filtered <- raw

# Rename columns with Sample IDs ------------------------------------------

names.info <- as.data.frame(names(raw.filtered)) %>% 
  select(sample=1) %>% 
  filter(str_detect(sample,'Geneid|Chr|Start|End|Strand|Length',negate = T)) %>% 
  mutate(bear = case_when(
    str_sub(sample,1,1) == 'C' ~ 'Cooke',
    str_sub(sample,1,1) == 'D' ~ 'Dodge',
    str_sub(sample,1,1) == 'F' ~ 'Frank',
    str_sub(sample,1,1) == 'Z' ~ 'Zuri',
    str_sub(sample,1,1) == 'J' ~ 'John',
    str_sub(sample,1,1) == 'O' ~ 'Oakley',
    str_sub(sample,1,1) == 'P' ~ 'Peeka',
    str_sub(sample,1,1) == 'R' ~ 'Roan'
  )) %>% 
  mutate(treatment = case_when(
    str_sub(sample,2,) == 'HH' ~ 'HH',
    str_sub(sample,2,) == 'HA' ~ 'HA',
    str_sub(sample,2,) == 'HPD' ~ 'HG',
    str_sub(sample,2,) == 'AA' ~ 'AA',
    str_sub(sample,2,) == 'AH' ~ 'AH',
    str_sub(sample,2,) == 'APD' ~ 'AG',
  ))  %>% 
  mutate(sex = case_when(
    bear %in% c('Kio','Cooke','Willow','Zuri','Luna','Oakley','Peeka') ~ 'Female',
    bear %in% c('John','Adak','Dodge','Frank','Roan') ~ 'Male'
  )) %>% 
  mutate(combined = paste(sample,bear,sex,treatment,sep = ':'))


# Split raw count table by tissue and add detailed sample names -----------

adipose.raw <- raw.filtered %>%
  pivot_longer(-c(1,2,3,4,5,6),names_to = 'orig_name',values_to = 'count') %>%
  left_join(names.info,by=c('orig_name'='sample')) %>%
  dplyr::select(Geneid,Length,combined,count) %>%
  pivot_wider(names_from = combined,values_from = count)



# Pre-filtering -----------------------------------------------------------
# Requiring non-zero count in at least 10 samples. Might revisit later.

adipose.raw.prefilt <- adipose.raw %>% filter(rowSums(.[,-c(1,2)] != 0 ) >= 10,)


# Convert to dataframe ----------------------------------------------------

adipose.raw.prefilt.df <- adipose.raw.prefilt %>% dplyr::select(-Length) %>% column_to_rownames('Geneid')


# Run DEseq2 - Adipose  ---------------------------------------------------

adip.treatment <- factor(str_split_fixed(names(adipose.raw.prefilt.df),'[:]',4)[,4])
adip.individual <- factor(str_split_fixed(names(adipose.raw.prefilt.df),'[:]',4)[,2])
adip.sex <- factor(str_split_fixed(names(adipose.raw.prefilt.df),'[:]',4)[,3])

adip.colData <- DataFrame(treatment = adip.treatment,individual = adip.individual,sex = adip.sex)

adip.dds <- DESeqDataSetFromMatrix(adipose.raw.prefilt.df,adip.colData,formula(~sex + treatment))

adip.dds <- DESeq(adip.dds)

adip.normcounts <- counts(adip.dds,normalized=TRUE)
adip.vsd.normCounts <- as.data.frame(assay(vst(adip.dds, blind=FALSE)))


### Vst PCA
adip.vsd <- vst(adip.dds, blind=FALSE)
# write_tsv(as.data.frame(assay(adip.vsd)) %>% rownames_to_column('gene_id'),'analysis/gene_level_rnaseq/norm_counts/cell_vstNormCounts_02.13.23.tsv')

adip.pcaData <- plotPCA(adip.vsd, intgroup=c("treatment","sex"), returnData=TRUE, ntop=5000)

adip.pcaData <- adip.pcaData %>%
  mutate(sample = str_split_fixed(name,'[:]',2)[,1])

adip.percentVar <- round(100 * attr(adip.pcaData, "percentVar"))

p.1 <- ggplot(subset(adip.pcaData,treatment %in% c('HG','HH')), aes(PC1, PC2, color=treatment, label=sample,shape=sex)) +
  geom_point(size=4) +
  ggrepel::geom_text_repel() +
  # ggrepel::geom_text_repel(data=subset(adip.pcaData,sample=='ZF1N')) + # Highlight potential sample swap
  xlab(paste0("PC1: ",adip.percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",adip.percentVar[2],"% variance")) +
  ggtitle("Adipose - Top 5000 Genes") +
  theme_linedraw(base_size = 16) + theme(plot.title.position = 'plot',plot.title = element_text(face = 'bold'))

p.1

###
### HPD vs. HH
###

adip.deRes.fed_vs_hib <- as.data.frame(results(adip.dds, contrast=c('treatment','HG','HH')))
adip.deRes.fed_vs_hib.ihwRes <- ihw(pvalue ~ baseMean,  data = adip.deRes.fed_vs_hib, alpha = 0.05)
rejections(adip.deRes.fed_vs_hib.ihwRes)
adip.deRes.fed_vs_hib$IHW_pvalue <- adip.deRes.fed_vs_hib.ihwRes@df$adj_pvalue
adip.deRes.fed_vs_hib <- adip.deRes.fed_vs_hib[order(adip.deRes.fed_vs_hib$IHW_pvalue),]
adip.deRes.fed_vs_hib.sig <- adip.deRes.fed_vs_hib %>% filter(IHW_pvalue < 0.05)%>%
  rownames_to_column('gene_id') %>% separate(gene_id,into = c('gene_number','id','symbol'),sep = '[:]')

adip.deRes.fed_vs_hib.sig %>% filter(log2FoldChange>0) %>% nrow()
adip.deRes.fed_vs_hib.sig %>% filter(log2FoldChange<0) %>% nrow()

# write_tsv(adip.deRes.fed_vs_hib.sig,'analysis/gene_level_rnaseq/pairwise_results/adipose_PWRes_Cells_HPD.vs.HH.tsv')


