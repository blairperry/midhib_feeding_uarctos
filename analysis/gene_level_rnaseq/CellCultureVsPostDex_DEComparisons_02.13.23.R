

library(tidyverse)
library(ggvenn)
library(patchwork)
library(org.Hs.eg.db)
# BiocManager::install("clusterProfiler")
library('clusterProfiler')
library(DOSE)

# Adipose -----------------------------------------------------------------

adi.deRes.fed_vs_notFed <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/adipose_PWRes_Fed.vs.NotFed.tsv')
adi.deRes.fed_vs_hib <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/adipose_PWRes_Fed.vs.PDHib.tsv')
cell.deRes.fed_vs_notFed <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/adipose_PWRes_Cells_HPD.vs.HH.tsv')



adi.venn <- list('Adipse - Post-Dex Fed vs.\nNot Fed' = adi.deRes.fed_vs_notFed$id,
                 'Cell - Hib Cell + PF Serum vs.\nHib Cell + Hib Serum' = cell.deRes.fed_vs_notFed$id)

p.adi2 <- ggvenn(adi.venn,set_name_size = 4,fill_color=c('white','white','white'),show_percentage = F) + ggtitle('Adipose') + theme(plot.title = element_text(face='bold',size=18))

p.adi2

adi.deRes.overlap <- adi.deRes.fed_vs_notFed %>%
  filter(gene_number %in% cell.deRes.fed_vs_notFed$gene_number) %>%
  left_join(cell.deRes.fed_vs_notFed,by=c('gene_number'),suffix = c('_Tissue_PostDex','_Cell_PostDex')) %>% 
  mutate(fill_color = ifelse((log2FoldChange_Tissue_PostDex < 0 & log2FoldChange_Cell_PostDex < 0) | (log2FoldChange_Tissue_PostDex > 0 & log2FoldChange_Cell_PostDex > 0),F,T))

p.adi <- ggplot(adi.deRes.overlap,aes(y=log2FoldChange_Tissue_PostDex,x=log2FoldChange_Cell_PostDex)) +
  geom_point(size=2,alpha=0.6,pch=21,fill='black',show.legend = F) +
  geom_hline(yintercept = 0,lwd=0.2) + geom_vline(xintercept = 0,lwd=0.2)+
  geom_abline(slope = 1,lty=2,alpha=0.5)+
  labs(title='',x='Cell - HH vs. HG Log2FC',y='Tissue - Not Fed vs. Fed Log2FC') +
  # coord_cartesian(xlim = c(-5,5),ylim=c(-3,3)) +
  theme_linedraw(base_size = 14) + theme(plot.title.position = 'plot',plot.title = element_text(face='bold'))
p.adi


adi.deRes.overlap.sameDir <- adi.deRes.overlap %>%
  filter((log2FoldChange_Tissue_PostDex * log2FoldChange_Cell_PostDex) > 0) # filter to genes with same directions
nrow(adi.deRes.overlap.sameDir) / nrow(adi.deRes.overlap) # ~68.5% of overlapping DE genes share direction


# Checking overlap with reversed genes ------------------------------------
adi.deRes.hib_vs_act <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/adipose_J19Hib.vs.J19Act.tsv')

adi.reverse <- adi.deRes.fed_vs_notFed %>% 
  filter(gene_number %in% adi.deRes.hib_vs_act$gene_number) %>% 
  left_join(adi.deRes.hib_vs_act,by=c('gene_number'),suffix = c('_PostDexExp','_Jansen2019')) %>% 
  filter((log2FoldChange_PostDexExp * log2FoldChange_Jansen2019) < 0) # filter to genes with opposite directions 

cell.deRes.fed_vs_notFed$id %in% adi.reverse$id_PostDexExp %>% sum()

cell.deRes.fed_vs_notFed.reversed <- cell.deRes.fed_vs_notFed %>% 
  filter(id %in% adi.reverse$id_PostDexExp)

adi.deRes.overlap.sameDir %>% filter(id_Tissue_PostDex %in% adi.reverse$id_PostDexExp)


