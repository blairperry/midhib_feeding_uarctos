
library(tidyverse)
library(colorspace)
library(pheatmap)

insr.candList <- read_tsv('analysis/insulin_genes/InsulinSignaling_Uniprot_CandGenes_AllSymbols.txt',col_names = 'gene')

tissue_colors <- c('Adipose'='#E69F00','Liver'='#56B4E9','Muscle'='#009E73')


# Intersect with ALL genes in transcriptome -------------------------------

gff.geneInfo <- read_tsv('data/0_reference/GCF_023065955.1_UrsArc1.0_genomic.gtf',col_names = F,comment = '#') %>% 
  filter(X3=='gene') %>% 
  dplyr::select(desc = 9) %>% 
  mutate(desc = str_split_fixed(desc,'[;]',2)[,1] %>% str_remove_all('ID=')) %>% 
  mutate(desc = str_remove_all(desc,'gene_id \"gene-|"')) %>% 
  mutate(present_in_bear = T) %>% 
  unique()

insr.candList.full <- read_tsv('analysis/insulin_genes/insulinReceptorSignalingPathway_Uniprot_10.20.21.tsv') %>% 
  janitor::clean_names() %>% 
  dplyr::select(entry_name,gene_names) %>% 
  mutate(entry_name = str_remove_all(entry_name,'_HUMAN')) %>% 
  mutate(gene_names = str_replace_all(gene_names,' ',';')) %>% 
  mutate(gene_names = paste(entry_name,';',gene_names,sep = '')) %>% 
  # mutate(gene_names = str_split(gene_names,';')) %>% 
  separate_rows(gene_names,sep = ';') %>% 
  dplyr::select(synonym = gene_names, main_gene = entry_name)
  
insr.candList.full.overlap <- insr.candList.full %>% 
  left_join(gff.geneInfo,by=c('synonym' = 'desc')) %>% 
  filter(present_in_bear == T) %>% 
  unique()

sum(insr.candList.full.overlap$main_gene %in% insr.candList.full$main_gene) / length(unique(insr.candList.full$main_gene)) 


# Intersect with season DE genes ------------------------------------------

adi.deRes.act_hib <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/adipose_J19Hib.vs.J19Act.tsv') %>% 
  filter(symbol %in% insr.candList.full.overlap$synonym) %>% mutate(tissue = 'adipose')
liv.deRes.act_hib <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/liver_J19Hib.vs.J19Act.tsv')%>% 
  filter(symbol %in% insr.candList.full.overlap$synonym)%>% mutate(tissue = 'liver')
mus.deRes.act_hib <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/muscle_J19Hib.vs.J19Act.tsv')%>% 
  filter(symbol %in% insr.candList.full.overlap$synonym)%>% mutate(tissue = 'muscle')

all.insr.deRes.act_hib <- adi.deRes.act_hib %>% 
  bind_rows(liv.deRes.act_hib,mus.deRes.act_hib)

ggplot(all.insr.deRes.act_hib,aes(x=log2FoldChange,y=symbol,fill=-log10(IHW_pvalue))) +
  geom_segment(aes(x=0,xend=log2FoldChange,y=symbol,yend=symbol)) +
  geom_point(pch=21,size=6) +
  facet_wrap(~tissue,ncol = 1,scales = 'free_y') +
  theme_linedraw()



# Intersect with gene-level DE --------------------------------------------

adi.deRes.fed_notFed <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/adipose_PWRes_Fed.vs.NotFed.tsv') %>% 
  filter(symbol %in% insr.candList.full.overlap$synonym) %>% mutate(tissue = 'adipose')
liv.deRes.fed_notFed <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/liver_PWRes_Fed.vs.NotFed.tsv')%>% 
  filter(symbol %in% insr.candList.full.overlap$synonym)%>% mutate(tissue = 'liver')
mus.deRes.fed_notFed <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/muscle_PWRes_Fed.vs.NotFed.tsv')%>% 
  filter(symbol %in% insr.candList.full.overlap$synonym)%>% mutate(tissue = 'muscle')

all.insr.deRes.fed_notFed <- adi.deRes.fed_notFed %>% 
  bind_rows(liv.deRes.fed_notFed,mus.deRes.fed_notFed)

ggplot(all.insr.deRes.fed_notFed,aes(x=log2FoldChange,y=symbol,fill=-log10(IHW_pvalue))) +
  geom_segment(aes(x=0,xend=log2FoldChange,y=symbol,yend=symbol)) +
  geom_point(pch=21,size=6) +
  facet_wrap(~tissue,ncol = 1,scales = 'free_y') +
  theme_linedraw()


insr.DEheatdata <- insr.candList.full.overlap %>% 
  mutate(adipose = ifelse(synonym %in% adi.deRes.fed_notFed$symbol,1,0)) %>% 
  mutate(liver = ifelse(synonym %in% liv.deRes.fed_notFed$symbol,1,0)) %>% 
  mutate(muscle = ifelse(synonym %in% mus.deRes.fed_notFed$symbol,1,0)) %>% 
  arrange(synonym) %>% 
  mutate(synonym = factor(synonym,levels = .$synonym))


pheatmap(as.data.frame(insr.DEheatdata %>% dplyr::select(synonym,adipose,liver,muscle)) %>% column_to_rownames('synonym'),
         cluster_rows = F, cluster_cols = F,border_color = 'black')


# Check if any are "reversed" genes

adi.revGenes <- read_tsv('analysis/gene_level_rnaseq/reversed_genes/adipose_reversedGenes_02.28.23.tsv')
liv.revGenes <- read_tsv('analysis/gene_level_rnaseq/reversed_genes/liver_reversedGenes_02.28.23.tsv')
mus.revGenes <- read_tsv('analysis/gene_level_rnaseq/reversed_genes/muscle_reversedGenes_02.28.23.tsv')

adi.deRes.fed_notFed %>% filter(symbol %in% adi.revGenes$symbol_PostDexExp) #NOT REVERSED
liv.deRes.fed_notFed %>% filter(symbol %in% liv.revGenes$symbol_PostDexExp) #NOT REVERSED
mus.deRes.fed_notFed %>% filter(symbol %in% mus.revGenes$symbol_PostDexExp) #NOT REVERSED




# Intersect w/ DIU results ------------------------------------------------

adi.diu.insr <- read_tsv('analysis/isoform_level_rnaseq/isoformSwitchResults/adipose_allSigSwitches_02.24.23.tsv') %>% 
  mutate(symbol = str_split_fixed(gene_id,'[-]',2)[,2]) %>% 
  mutate(insrGene = ifelse(symbol %in% insr.candList$gene,T,F)) %>% 
  mutate(tissue = 'Adipose',result='Differential isoform usage')

liv.diu.insr <- read_tsv('analysis/isoform_level_rnaseq/isoformSwitchResults/liver_allSigSwitches_02.24.23.tsv') %>% 
  mutate(symbol = str_split_fixed(gene_id,'[-]',2)[,2]) %>% 
  mutate(insrGene = ifelse(symbol %in% insr.candList$gene,T,F)) %>% 
  mutate(tissue = 'Liver',result='Differential isoform usage')

mus.diu.insr <- read_tsv('analysis/isoform_level_rnaseq/isoformSwitchResults/muscle_allSigSwitches_02.24.23.tsv') %>% 
  mutate(symbol = str_split_fixed(gene_id,'[-]',2)[,2]) %>% 
  mutate(insrGene = ifelse(symbol %in% insr.candList$gene,T,F)) %>% 
  mutate(tissue = 'Muscle',result='Differential isoform usage')

all.diu.insr <- adi.diu.insr %>% 
  bind_rows(liv.diu.insr,mus.diu.insr) %>% 
  filter(insrGene==T)

