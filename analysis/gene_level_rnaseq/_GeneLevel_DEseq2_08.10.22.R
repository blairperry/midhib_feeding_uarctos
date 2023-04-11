
library(tidyverse)
library(DESeq2)
library(IHW)
library(patchwork)

# Read in and reformat raw count table ------------------------------------

raw <- read_tsv('data/3_featureCounts/postDexExperiment_08.10.22.txt',comment = '#')

# Clean sample names
names(raw) <- str_remove_all(names(raw),'/data/kelley/projects/bear_dextrose_reversal_newAnalysis_bwp/2_star/star_mapped/|Aligned.sortedByCoord.out.bam')

# Remove samples with fewer than 3M reads mapped, known contaminants, known outlier samples, and CFA (previously shown to contain hair follicle)
raw.filtered <- raw %>%
  dplyr::select(-CFA_S9,-D2_S26,-H2_S57,-E2_S34,-D7_S31,-H7_S62) %>% # remove CFA and low mapping
  dplyr::select(-C2_S18,-E3_S35,-F2_S42) %>%  # remove contaminated samples
  dplyr::select(-F3_S43,-H6_S61,-A2_S2) # outlier samples

# Rename columns with Sample IDs ------------------------------------------

names.temp <- as.data.frame(names(raw.filtered)) %>%
  dplyr::select(orig_name = 1) %>%
  filter(str_detect(orig_name,'[_]')) %>%
  mutate(simple = str_split_fixed(orig_name,'[_]',2)[,1]) %>%
  mutate(simple = ifelse(str_detect(simple,'Unknown'),orig_name,simple))

sample.info <- readxl::read_xlsx('data/sampleData.xlsx') %>% janitor::clean_names() %>% dplyr::select(1:3,season=4)
sample.info.more <- readxl::read_xlsx('data/BearPostDexLibrary Information 2019.xlsx') %>% janitor::clean_names() %>% dplyr::select(1,4:7)

names.info <- names.temp %>%
  left_join(sample.info,by=c('simple'='barcode')) %>%
  dplyr::select(-tissue) %>%
  left_join(sample.info.more,by=c('sample'='sample_tube')) %>%
  mutate(sample = ifelse(is.na(sample),simple,sample)) %>%
  mutate(tissue = ifelse(is.na(tissue),
                         case_when(
                           str_sub(sample,2,2) == 'F' ~ 'Fat',
                           str_sub(sample,2,2) == 'L' ~ 'Liver',
                           str_sub(sample,2,2) == 'M' ~ 'Muscle'
                         ),
                         tissue)) %>%
  mutate(tissue = ifelse(is.na(tissue),'Fat',tissue)) %>%
  mutate(bear = ifelse(is.na(bear),
                         case_when(
                           str_sub(sample,1,1) == 'C' ~ 'Cooke',
                           str_sub(sample,1,1) == 'D' ~ 'Dodge',
                           str_sub(sample,1,1) == 'F' ~ 'Frank',
                           str_sub(sample,1,1) == 'Z' ~ 'Zuri',
                           str_sub(sample,1,1) == 'J' ~ 'John',
                           str_sub(sample,1,1) == 'O' ~ 'Oakley',
                           str_sub(sample,1,1) == 'P' ~ 'Peeka',
                           str_sub(sample,1,1) == 'R' ~ 'Roan'
                         ),
                       bear)) %>%
  mutate(bear = ifelse(is.na(bear),sample,bear)) %>%
  mutate(treatment = case_when(
    x1st_or_2nd_sampling_period == '1st' ~ 'PDExp_Hib',
    x1st_or_2nd_sampling_period == '2nd'& fed_dextrose == 'Yes' ~ 'PDExp_Fed',
    x1st_or_2nd_sampling_period == '2nd'& fed_dextrose == 'No' ~ 'PDExp_NotFed',
    str_detect(sample,'Hi') ~ 'Jansen2019_Hib',
    str_detect(sample,'Hy') ~ 'Jansen2019_Hyp',
    sample %in% c("Adak","Cooke","Dodge","Frank","John","Kio","Oakley","Peeka","Roan","Unknown_CAGATC","Unknown_TGACCA","Zuri") ~ 'PDExp_50Per'
    )) %>%
  mutate(treatment = ifelse(is.na(treatment),'Jansen2019_Act',treatment)) %>%
  mutate(sex = case_when(
    bear %in% c('Kio','Cooke','Willow','Zuri','Luna','Oakley','Peeka') ~ 'Female',
    bear %in% c('John','Adak','Dodge','Frank','Roan') ~ 'Male',
    str_detect(bear,'Unkn') ~ 'Unknown'
  )) %>%
  dplyr::select(orig_name,sample,tissue,bear,sex,treatment) %>%
  mutate(combined = paste(sample,tissue,bear,sex,treatment,sep = ':'))

# write_tsv(names.info,'sampleInfo_08.11.22.txt')

# Split raw count table by tissue and add detailed sample names -----------

adipose.raw <- raw.filtered %>%
  pivot_longer(-c(1,2,3,4,5,6),names_to = 'orig_name',values_to = 'count') %>%
  left_join(names.info) %>%
  filter(tissue == 'Fat') %>%
  filter(str_detect(bear,'Unknown',negate = T)) %>%  # Filtering out unknown 50% samples for now, since we don't know sex
  dplyr::select(Geneid,Length,combined,count) %>%
  pivot_wider(names_from = combined,values_from = count)

liver.raw <- raw.filtered %>%
  pivot_longer(-c(1,2,3,4,5,6),names_to = 'orig_name',values_to = 'count') %>%
  left_join(names.info) %>%
  filter(tissue == 'Liver') %>%
  filter(str_detect(bear,'Unknown',negate = T)) %>%  # Filtering out unknown 50% samples for now, since we don't know sex
  dplyr::select(Geneid,Length,combined,count) %>%
  pivot_wider(names_from = combined,values_from = count)

muscle.raw <- raw.filtered %>%
  pivot_longer(-c(1,2,3,4,5,6),names_to = 'orig_name',values_to = 'count') %>%
  left_join(names.info) %>%
  filter(tissue == 'Muscle') %>%
  filter(str_detect(bear,'Unknown',negate = T)) %>%  # Filtering out unknown 50% samples for now, since we don't know sex
  dplyr::select(Geneid,Length,combined,count) %>%
  pivot_wider(names_from = combined,values_from = count)


# Pre-filtering -----------------------------------------------------------
# Requiring non-zero count in at least 10 samples.

adipose.raw.prefilt <- adipose.raw %>% filter(rowSums(.[,-c(1,2)] != 0 ) >= 10,)
liver.raw.prefilt <- liver.raw %>% filter(rowSums(.[,-c(1,2)] != 0 ) >= 10,)
muscle.raw.prefilt <- muscle.raw %>% filter(rowSums(.[,-c(1,2)] != 0 ) >= 10,)


# Convert to dataframe ----------------------------------------------------

adipose.raw.prefilt.df <- adipose.raw.prefilt %>% dplyr::select(-Length) %>% column_to_rownames('Geneid')
liver.raw.prefilt.df <- liver.raw.prefilt  %>% dplyr::select(-Length) %>% column_to_rownames('Geneid')
muscle.raw.prefilt.df <- muscle.raw.prefilt %>% dplyr::select(-Length) %>% column_to_rownames('Geneid')


# Run DEseq2 - Adipose  ---------------------------------------------------

adip.treatment <- factor(str_split_fixed(names(adipose.raw.prefilt.df),'[:]',5)[,5])
adip.individual <- factor(str_split_fixed(names(adipose.raw.prefilt.df),'[:]',5)[,3])
adip.sex <- factor(str_split_fixed(names(adipose.raw.prefilt.df),'[:]',5)[,4])

adip.colData <- DataFrame(treatment = adip.treatment,individual = adip.individual,sex = adip.sex)

adip.dds <- DESeqDataSetFromMatrix(adipose.raw.prefilt.df,adip.colData,formula(~sex + treatment))

adip.dds <- DESeq(adip.dds)

adip.normcounts <- counts(adip.dds,normalized=TRUE)
adip.vsd.normCounts <- as.data.frame(assay(vst(adip.dds, blind=FALSE)))

### Vst PCA
adip.vsd <- vst(adip.dds, blind=FALSE)
# write_tsv(as.data.frame(assay(adip.vsd)) %>% rownames_to_column('gene_id'),'analysis/gene_level_rnaseq/norm_counts/adipose_vstNormCounts_08.11.22.tsv')

adip.pcaData <- plotPCA(adip.vsd, intgroup=c("treatment","sex"), returnData=TRUE, ntop=5000)

adip.pcaData <- adip.pcaData %>%
  mutate(sample = str_split_fixed(name,'[:]',2)[,1])

adip.percentVar <- round(100 * attr(adip.pcaData, "percentVar"))

p.1 <- ggplot(subset(adip.pcaData,treatment %in% c('Jansen2019_Act','PDExp_Hib','PDExp_Fed')), aes(PC1, PC2, color=treatment, label=sample,shape=sex)) +
  geom_point(size=4) +
  ggrepel::geom_text_repel() +
  # ggrepel::geom_text_repel(data=subset(adip.pcaData,sample=='ZF1N')) + # Highlight potential sample swap
  xlab(paste0("PC1: ",adip.percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",adip.percentVar[2],"% variance")) +
  ggtitle("Adipose - Top 5000 Genes") +
  theme_linedraw(base_size = 16) + theme(plot.title.position = 'plot',plot.title = element_text(face = 'bold'))

p.1


###
### PDExp_Fed vs. PDExp_NotFed
###

adip.deRes.fed_vs_notFed <- as.data.frame(results(adip.dds, contrast=c('treatment','PDExp_Fed','PDExp_NotFed')))
adip.deRes.fed_vs_notFed.ihwRes <- ihw(pvalue ~ baseMean,  data = adip.deRes.fed_vs_notFed, alpha = 0.05)
rejections(adip.deRes.fed_vs_notFed.ihwRes)
adip.deRes.fed_vs_notFed$IHW_pvalue <- adip.deRes.fed_vs_notFed.ihwRes@df$adj_pvalue
adip.deRes.fed_vs_notFed <- adip.deRes.fed_vs_notFed[order(adip.deRes.fed_vs_notFed$IHW_pvalue),]
adip.deRes.fed_vs_notFed.sig <- adip.deRes.fed_vs_notFed %>% filter(IHW_pvalue < 0.05)%>%
  rownames_to_column('gene_id') %>% separate(gene_id,into = c('gene_number','id','symbol'),sep = '[:]')

nrow(adip.deRes.fed_vs_notFed.sig)
nrow(adip.deRes.fed_vs_notFed.sig %>% filter(log2FoldChange>0))
nrow(adip.deRes.fed_vs_notFed.sig %>% filter(log2FoldChange<0))

# write_tsv(adip.deRes.fed_vs_notFed.sig,'analysis/gene_level_rnaseq/pairwise_results/adipose_PWRes_Fed.vs.NotFed.tsv')

###
### PDExp_Fed vs. PDExp_Hib
###

adip.deRes.fed_vs_pdhib <- as.data.frame(results(adip.dds, contrast=c('treatment','PDExp_Fed','PDExp_Hib')))
adip.deRes.fed_vs_pdhib.ihwRes <- ihw(pvalue ~ baseMean,  data = adip.deRes.fed_vs_pdhib, alpha = 0.05)
rejections(adip.deRes.fed_vs_pdhib.ihwRes)
adip.deRes.fed_vs_pdhib$IHW_pvalue <- adip.deRes.fed_vs_pdhib.ihwRes@df$adj_pvalue
adip.deRes.fed_vs_pdhib <- adip.deRes.fed_vs_pdhib[order(adip.deRes.fed_vs_pdhib$IHW_pvalue),]
adip.deRes.fed_vs_pdhib.sig <- adip.deRes.fed_vs_pdhib %>% filter(IHW_pvalue < 0.05)%>%
  rownames_to_column('gene_id') %>% separate(gene_id,into = c('gene_number','id','symbol'),sep = '[:]')

nrow(adip.deRes.fed_vs_pdhib.sig)
nrow(adip.deRes.fed_vs_pdhib.sig %>% filter(log2FoldChange>0))
nrow(adip.deRes.fed_vs_pdhib.sig %>% filter(log2FoldChange<0))


# write_tsv(adip.deRes.fed_vs_pdhib.sig,'analysis/gene_level_rnaseq/pairwise_results/adipose_PWRes_Fed.vs.PDHib.tsv')

###
### PDExp_NotFed vs. PDExp_Hib
###

adip.deRes.notFed_vs_pdhib <- as.data.frame(results(adip.dds, contrast=c('treatment','PDExp_NotFed','PDExp_Hib')))
adip.deRes.fed_vs_notFed.ihwRes <- ihw(pvalue ~ baseMean,  data = adip.deRes.notFed_vs_pdhib, alpha = 0.05)
rejections(adip.deRes.fed_vs_notFed.ihwRes) # None significant

###
### Jansen2019_Hib vs. Jansen2019_Act
###

adip.deRes.hib_vs_act <- as.data.frame(results(adip.dds, contrast=c('treatment','Jansen2019_Hib','Jansen2019_Act')))
adip.deRes.hib_vs_act.ihwRes <- ihw(pvalue ~ baseMean,  data = adip.deRes.hib_vs_act, alpha = 0.05)
rejections(adip.deRes.hib_vs_act.ihwRes)
adip.deRes.hib_vs_act$IHW_pvalue <- adip.deRes.hib_vs_act.ihwRes@df$adj_pvalue
adip.deRes.hib_vs_act <- adip.deRes.hib_vs_act[order(adip.deRes.hib_vs_act$IHW_pvalue),]
adip.deRes.hib_vs_act.sig <- adip.deRes.hib_vs_act %>% filter(IHW_pvalue < 0.05)%>%
  rownames_to_column('gene_id') %>% separate(gene_id,into = c('gene_number','id','symbol'),sep = '[:]')

nrow(adip.deRes.hib_vs_act.sig)
nrow(adip.deRes.hib_vs_act.sig %>% filter(log2FoldChange>0))
nrow(adip.deRes.hib_vs_act.sig %>% filter(log2FoldChange<0))


# write_tsv(adip.deRes.hib_vs_act.sig,'analysis/gene_level_rnaseq/pairwise_results/adipose_J19Hib.vs.J19Act.tsv')



# Run DEseq2 - Liver  ---------------------------------------------------

liv.treatment <- factor(str_split_fixed(names(liver.raw.prefilt.df),'[:]',5)[,5])
liv.individual <- factor(str_split_fixed(names(liver.raw.prefilt.df),'[:]',5)[,3])
liv.sex <- factor(str_split_fixed(names(liver.raw.prefilt.df),'[:]',5)[,4])

liv.colData <- DataFrame(treatment = liv.treatment,individual = liv.individual,sex = liv.sex)

liv.dds <- DESeqDataSetFromMatrix(liver.raw.prefilt.df,liv.colData,formula(~sex + treatment))

liv.dds <- DESeq(liv.dds)

liv.normcounts <- counts(liv.dds,normalized=TRUE)
liv.vsd.normCounts <- as.data.frame(assay(vst(liv.dds, blind=FALSE)))

### Vst PCA
liv.vsd <- vst(liv.dds, blind=FALSE)
# write_tsv(as.data.frame(assay(liv.vsd)) %>% rownames_to_column('gene_id'),'analysis/gene_level_rnaseq/norm_counts/liver_vstNormCounts_08.11.22.tsv')

liv.pcaData <- plotPCA(liv.vsd, intgroup=c("treatment","sex"), returnData=TRUE, ntop=5000)

liv.pcaData <- liv.pcaData %>%
  mutate(sample = str_split_fixed(name,'[:]',2)[,1])

liv.percentVar <- round(100 * attr(liv.pcaData, "percentVar"))

p.2 <- ggplot(subset(liv.pcaData,treatment %in% c('Jansen2019_Act','PDExp_Hib','PDExp_Fed')), aes(PC1, PC2, color=treatment, label=sample,shape=sex)) +
  geom_point(size=4) +
  # ggrepel::geom_text_repel() +
  xlab(paste0("PC1: ",liv.percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",liv.percentVar[2],"% variance")) +
  ggtitle("Liver - Top 5000 Genes") +
  theme_linedraw(base_size = 16) + theme(plot.title.position = 'plot',plot.title = element_text(face = 'bold'))


###
### PDExp_Fed vs. PDExp_NotFed
###

liv.deRes.fed_vs_notFed <- as.data.frame(results(liv.dds, contrast=c('treatment','PDExp_Fed','PDExp_NotFed')))
liv.deRes.fed_vs_notFed.ihwRes <- ihw(pvalue ~ baseMean,  data = liv.deRes.fed_vs_notFed, alpha = 0.05)
rejections(liv.deRes.fed_vs_notFed.ihwRes)
liv.deRes.fed_vs_notFed$IHW_pvalue <- liv.deRes.fed_vs_notFed.ihwRes@df$adj_pvalue
liv.deRes.fed_vs_notFed <- liv.deRes.fed_vs_notFed[order(liv.deRes.fed_vs_notFed$IHW_pvalue),]
liv.deRes.fed_vs_notFed.sig <- liv.deRes.fed_vs_notFed %>% filter(IHW_pvalue < 0.05)%>%
  rownames_to_column('gene_id') %>% separate(gene_id,into = c('gene_number','id','symbol'),sep = '[:]')

nrow(liv.deRes.fed_vs_notFed.sig)
nrow(liv.deRes.fed_vs_notFed.sig %>% filter(log2FoldChange>0))
nrow(liv.deRes.fed_vs_notFed.sig %>% filter(log2FoldChange<0))


liv.deRes.fed_vs_notFed %>% rownames_to_column() %>% filter(str_detect(rowname,'FOXO1'))

# write_tsv(liv.deRes.fed_vs_notFed.sig,'analysis/gene_level_rnaseq/pairwise_results/liver_PWRes_Fed.vs.NotFed.tsv')


###
### PDExp_Fed vs. PDExp_Hib
###

liv.deRes.fed_vs_pdhib <- as.data.frame(results(liv.dds, contrast=c('treatment','PDExp_Fed','PDExp_Hib')))
liv.deRes.fed_vs_pdhib.ihwRes <- ihw(pvalue ~ baseMean,  data = liv.deRes.fed_vs_pdhib, alpha = 0.05)
rejections(liv.deRes.fed_vs_pdhib.ihwRes)
liv.deRes.fed_vs_pdhib$IHW_pvalue <- liv.deRes.fed_vs_pdhib.ihwRes@df$adj_pvalue
liv.deRes.fed_vs_pdhib <- liv.deRes.fed_vs_pdhib[order(liv.deRes.fed_vs_pdhib$IHW_pvalue),]
liv.deRes.fed_vs_pdhib.sig <- liv.deRes.fed_vs_pdhib %>% filter(IHW_pvalue < 0.05)%>%
  rownames_to_column('gene_id') %>% separate(gene_id,into = c('gene_number','id','symbol'),sep = '[:]')

nrow(liv.deRes.fed_vs_pdhib.sig)
nrow(liv.deRes.fed_vs_pdhib.sig %>% filter(log2FoldChange>0))
nrow(liv.deRes.fed_vs_pdhib.sig %>% filter(log2FoldChange<0))

# write_tsv(liv.deRes.fed_vs_pdhib.sig,'analysis/gene_level_rnaseq/pairwise_results/liver_PWRes_Fed.vs.PDHib.tsv')

###
### PDExp_NotFed vs. PDExp_Hib
###

liv.deRes.notFed_vs_pdhib <- as.data.frame(results(liv.dds, contrast=c('treatment','PDExp_NotFed','PDExp_Hib')))
liv.deRes.notFed_vs_pdhib.ihwRes <- ihw(pvalue ~ baseMean,  data = liv.deRes.notFed_vs_pdhib, alpha = 0.05)
rejections(liv.deRes.notFed_vs_pdhib.ihwRes) # None significant
liv.deRes.notFed_vs_pdhib$IHW_pvalue <- liv.deRes.notFed_vs_pdhib.ihwRes@df$adj_pvalue
liv.deRes.notFed_vs_pdhib <- liv.deRes.notFed_vs_pdhib[order(liv.deRes.notFed_vs_pdhib$IHW_pvalue),]
liv.deRes.notFed_vs_pdhib.sig <- liv.deRes.notFed_vs_pdhib %>% filter(IHW_pvalue < 0.05)%>%
  rownames_to_column('gene_id') %>% separate(gene_id,into = c('gene_number','id','symbol'),sep = '[:]')

nrow(liv.deRes.notFed_vs_pdhib.sig)
nrow(liv.deRes.notFed_vs_pdhib.sig %>% filter(log2FoldChange>0))
nrow(liv.deRes.notFed_vs_pdhib.sig %>% filter(log2FoldChange<0))

# write_tsv(liv.deRes.notFed_vs_pdhib.sig,'analysis/gene_level_rnaseq/pairwise_results/liver_PWRes_NotFed.vs.PDHib.tsv')

###
### Jansen2019_Hib vs. Jansen2019_Act
###

liv.deRes.hib_vs_act <- as.data.frame(results(liv.dds, contrast=c('treatment','Jansen2019_Hib','Jansen2019_Act')))
liv.deRes.hib_vs_act.ihwRes <- ihw(pvalue ~ baseMean,  data = liv.deRes.hib_vs_act, alpha = 0.05)
rejections(liv.deRes.hib_vs_act.ihwRes)
liv.deRes.hib_vs_act$IHW_pvalue <- liv.deRes.hib_vs_act.ihwRes@df$adj_pvalue
liv.deRes.hib_vs_act <- liv.deRes.hib_vs_act[order(liv.deRes.hib_vs_act$IHW_pvalue),]
liv.deRes.hib_vs_act.sig <- liv.deRes.hib_vs_act %>% filter(IHW_pvalue < 0.05)%>%
  rownames_to_column('gene_id') %>% separate(gene_id,into = c('gene_number','id','symbol'),sep = '[:]')

nrow(liv.deRes.hib_vs_act.sig)
nrow(liv.deRes.hib_vs_act.sig %>% filter(log2FoldChange>0))
nrow(liv.deRes.hib_vs_act.sig %>% filter(log2FoldChange<0))

# write_tsv(liv.deRes.hib_vs_act.sig,'analysis/gene_level_rnaseq/pairwise_results/liver_J19Hib.vs.J19Act.tsv')




# Run DEseq2 - muscle  ---------------------------------------------------

mus.treatment <- factor(str_split_fixed(names(muscle.raw.prefilt.df),'[:]',5)[,5])
mus.individual <- factor(str_split_fixed(names(muscle.raw.prefilt.df),'[:]',5)[,3])
mus.sex <- factor(str_split_fixed(names(muscle.raw.prefilt.df),'[:]',5)[,4])

mus.colData <- DataFrame(treatment = mus.treatment,individual = mus.individual,sex = mus.sex)

mus.dds <- DESeqDataSetFromMatrix(muscle.raw.prefilt.df,mus.colData,formula(~sex + treatment))

mus.dds <- DESeq(mus.dds)

mus.normcounts <- counts(mus.dds,normalized=TRUE)
mus.vsd.normCounts <- as.data.frame(assay(vst(mus.dds, blind=FALSE)))

### Vst PCA
mus.vsd <- vst(mus.dds, blind=FALSE)
# write_tsv(as.data.frame(assay(mus.vsd)) %>% rownames_to_column('gene_id'),'analysis/gene_level_rnaseq/norm_counts/muscle_vstNormCounts_08.11.22.tsv')

mus.pcaData <- plotPCA(mus.vsd, intgroup=c("treatment","sex"), returnData=TRUE, ntop=5000)

mus.pcaData <- mus.pcaData %>%
  mutate(sample = str_split_fixed(name,'[:]',2)[,1]) %>%
  mutate(plot_treatment = case_when(
    str_detect(treatment,'Hib') ~ 'Hibernation',
    str_detect(treatment,'Act') ~ 'Active',
    str_detect(treatment,'NotFed') ~ 'Hibernation',
    str_detect(treatment,'Exp_Fed') ~ 'Fed during hibernation'
  )) %>%
  filter(!is.na(plot_treatment))


mus.percentVar <- round(100 * attr(mus.pcaData, "percentVar"))

p.3 <- ggplot(subset(mus.pcaData,treatment %in% c('Jansen2019_Act','PDExp_Hib','PDExp_Fed')), aes(PC1, PC2, color=treatment, label=sample,shape=sex)) +
  geom_point(size=4) +
  # ggrepel::geom_text_repel() +
  xlab(paste0("PC1: ",mus.percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",mus.percentVar[2],"% variance")) +
  ggtitle("Muscle - Top 5000 Genes") +
  theme_linedraw(base_size = 16) + theme(plot.title.position = 'plot',plot.title = element_text(face = 'bold'))


###
### PDExp_Fed vs. PDExp_NotFed
###

mus.deRes.fed_vs_notFed <- as.data.frame(results(mus.dds, contrast=c('treatment','PDExp_Fed','PDExp_NotFed')))
mus.deRes.fed_vs_notFed.ihwRes <- ihw(pvalue ~ baseMean,  data = mus.deRes.fed_vs_notFed, alpha = 0.05)
rejections(mus.deRes.fed_vs_notFed.ihwRes)
mus.deRes.fed_vs_notFed$IHW_pvalue <- mus.deRes.fed_vs_notFed.ihwRes@df$adj_pvalue
mus.deRes.fed_vs_notFed <- mus.deRes.fed_vs_notFed[order(mus.deRes.fed_vs_notFed$IHW_pvalue),]
mus.deRes.fed_vs_notFed.sig <- mus.deRes.fed_vs_notFed %>% filter(IHW_pvalue < 0.05)%>%
  rownames_to_column('gene_id') %>% separate(gene_id,into = c('gene_number','id','symbol'),sep = '[:]')

nrow(mus.deRes.fed_vs_notFed.sig)
nrow(mus.deRes.fed_vs_notFed.sig %>% filter(log2FoldChange>0))
nrow(mus.deRes.fed_vs_notFed.sig %>% filter(log2FoldChange<0))

# write_tsv(mus.deRes.fed_vs_notFed.sig,'analysis/gene_level_rnaseq/pairwise_results/muscle_PWRes_Fed.vs.NotFed.tsv')

###
### PDExp_Fed vs. PDExp_Hib
###

mus.deRes.fed_vs_pdhib <- as.data.frame(results(mus.dds, contrast=c('treatment','PDExp_Fed','PDExp_Hib')))
mus.deRes.fed_vs_pdhib.ihwRes <- ihw(pvalue ~ baseMean,  data = mus.deRes.fed_vs_pdhib, alpha = 0.05)
rejections(mus.deRes.fed_vs_pdhib.ihwRes)
mus.deRes.fed_vs_pdhib$IHW_pvalue <- mus.deRes.fed_vs_pdhib.ihwRes@df$adj_pvalue
mus.deRes.fed_vs_pdhib <- mus.deRes.fed_vs_pdhib[order(mus.deRes.fed_vs_pdhib$IHW_pvalue),]
mus.deRes.fed_vs_pdhib.sig <- mus.deRes.fed_vs_pdhib %>% filter(IHW_pvalue < 0.05) %>%
  rownames_to_column('gene_id') %>% separate(gene_id,into = c('gene_number','id','symbol'),sep = '[:]')

nrow(mus.deRes.fed_vs_pdhib.sig)
nrow(mus.deRes.fed_vs_pdhib.sig %>% filter(log2FoldChange>0))
nrow(mus.deRes.fed_vs_pdhib.sig %>% filter(log2FoldChange<0))

# write_tsv(mus.deRes.fed_vs_pdhib.sig,'analysis/gene_level_rnaseq/pairwise_results/muscle_PWRes_Fed.vs.PDHib.tsv')

###
### PDExp_NotFed vs. PDExp_Hib
###

mus.deRes.notFed_vs_pdhib <- as.data.frame(results(mus.dds, contrast=c('treatment','PDExp_NotFed','PDExp_Hib')))
mus.deRes.notFed_vs_pdhib.ihwRes <- ihw(pvalue ~ baseMean,  data = mus.deRes.notFed_vs_pdhib, alpha = 0.05)
rejections(mus.deRes.notFed_vs_pdhib.ihwRes) # None significant

###
### Jansen2019_Hib vs. Jansen2019_Act
###

mus.deRes.hib_vs_act <- as.data.frame(results(mus.dds, contrast=c('treatment','Jansen2019_Hib','Jansen2019_Act')))
mus.deRes.hib_vs_act.ihwRes <- ihw(pvalue ~ baseMean,  data = mus.deRes.hib_vs_act, alpha = 0.05)
rejections(mus.deRes.hib_vs_act.ihwRes)
mus.deRes.hib_vs_act$IHW_pvalue <- mus.deRes.hib_vs_act.ihwRes@df$adj_pvalue
mus.deRes.hib_vs_act <- mus.deRes.hib_vs_act[order(mus.deRes.hib_vs_act$IHW_pvalue),]
mus.deRes.hib_vs_act.sig <- mus.deRes.hib_vs_act %>% filter(IHW_pvalue < 0.05)%>%
  rownames_to_column('gene_id') %>% separate(gene_id,into = c('gene_number','id','symbol'),sep = '[:]')

nrow(mus.deRes.hib_vs_act.sig)
nrow(mus.deRes.hib_vs_act.sig %>% filter(log2FoldChange>0))
nrow(mus.deRes.hib_vs_act.sig %>% filter(log2FoldChange<0))


# write_tsv(mus.deRes.hib_vs_act.sig,'analysis/gene_level_rnaseq/pairwise_results/muscle_J19Hib.vs.J19Act.tsv')



