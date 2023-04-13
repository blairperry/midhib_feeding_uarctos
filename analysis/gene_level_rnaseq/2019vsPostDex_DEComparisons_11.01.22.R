
library(tidyverse)
library(ggvenn)
library(patchwork)
library(org.Hs.eg.db)
# BiocManager::install("clusterProfiler")
library('clusterProfiler')
library(DOSE)



# Adipose -----------------------------------------------------------------

adi.deRes.fed_vs_notFed <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/adipose_PWRes_Fed.vs.NotFed.tsv')
adi.deRes.fed_vs_pdhib <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/adipose_PWRes_Fed.vs.PDHib.tsv')
adi.deRes.hib_vs_act <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/adipose_J19Hib.vs.J19Act.tsv')

adi.venn <- list('Post-Dex Fed vs.\nPost-Dex Not Fed' = adi.deRes.fed_vs_notFed$id,
                 # 'Post-Dex Fed vs.\nPost-Dex Hib' = adi.deRes.fed_vs_pdhib$id,
                 'Jansen 2019 Hib vs.\nJansen 2019 Act' = adi.deRes.hib_vs_act$id)


p.adi2 <- ggvenn(adi.venn,set_name_size = 4,fill_color=c('white','white','white'),show_percentage = F) + ggtitle('Adipose') + theme(plot.title = element_text(face='bold',size=18))

adi.deRes.overlap <- adi.deRes.fed_vs_notFed %>% 
  filter(gene_number %in% adi.deRes.hib_vs_act$gene_number) %>% 
  left_join(adi.deRes.hib_vs_act,by=c('gene_number'),suffix = c('_PostDexExp','_Jansen2019'))

p.adi <- ggplot(adi.deRes.overlap,aes(y=log2FoldChange_PostDexExp,x=log2FoldChange_Jansen2019)) +
  geom_point(size=2,alpha=0.6,col='dodgerblue4') +
  geom_hline(yintercept = 0,lwd=0.2) + geom_vline(xintercept = 0,lwd=0.2)+
  geom_abline(slope = -1,lty=2,alpha=0.5)+  
  labs(title='',x='Active vs. Hibernation Log2FC',y='Not Fed vs. Fed Log2FC') +
  # coord_cartesian(xlim = c(-5,5),ylim=c(-3,3)) +
  theme_linedraw(base_size = 14) + theme(plot.title.position = 'plot',plot.title = element_text(face='bold'))

adi.deRes.overlap.oppDir <- adi.deRes.overlap %>% 
  filter((log2FoldChange_PostDexExp * log2FoldChange_Jansen2019) < 0) # filter to genes with opposite directions 
# write_tsv(adi.deRes.overlap.oppDir,'analysis/gene_level_rnaseq/reversed_genes/adipose_reversedGenes_02.28.23.tsv')

adi.deRes.overlap.oppDir %>% filter(log2FoldChange_PostDexExp > 0) %>% nrow()

adi.deRes.pdexpUnique <- adi.deRes.fed_vs_notFed %>% 
  filter(!(gene_number %in% adi.deRes.hib_vs_act$gene_number)) 
adi.deRes.jansenUnique <- adi.deRes.hib_vs_act %>% 
  filter(!(gene_number %in% adi.deRes.fed_vs_notFed$gene_number)) 


# Liver -------------------------------------------------------------------

liv.deRes.fed_vs_notFed <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/liver_PWRes_Fed.vs.NotFed.tsv')
liv.deRes.fed_vs_pdhib <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/liver_PWRes_Fed.vs.PDHib.tsv')
liv.deRes.hib_vs_act <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/liver_J19Hib.vs.J19Act.tsv')

liv.venn <- list('Post-Dex Fed vs.\nPost-Dex Not Fed' = liv.deRes.fed_vs_notFed$id,
                 # 'Post-Dex Fed vs.\nPost-Dex Hib' = liv.deRes.fed_vs_pdhib$id,
                 'Jansen 2019 Hib vs.\nJansen 2019 Act' = liv.deRes.hib_vs_act$id)

p.liv2 <- ggvenn(liv.venn,set_name_size = 4,fill_color=c('white','white','white'),show_percentage = F) + ggtitle('Liver') + theme(plot.title = element_text(face='bold',size=18))

liv.deRes.overlap <- liv.deRes.fed_vs_notFed %>% 
  filter(gene_number %in% liv.deRes.hib_vs_act$gene_number) %>% 
  left_join(liv.deRes.hib_vs_act,by=c('gene_number'),suffix = c('_PostDexExp','_Jansen2019'))

p.liv <- ggplot(liv.deRes.overlap,aes(y=log2FoldChange_PostDexExp,x=log2FoldChange_Jansen2019)) +
  geom_point(size=2,alpha=0.6,col='goldenrod') +
  geom_hline(yintercept = 0,lwd=0.2) + geom_vline(xintercept = 0,lwd=0.2)+
  geom_abline(slope = -1,lty=2,alpha=0.5)+  
  labs(title='',x='Active vs. Hibernation Log2FC',y=' ') +
  # coord_cartesian(xlim = c(-5,5),ylim=c(-3,3)) +
  theme_linedraw(base_size = 14) + theme(plot.title.position = 'plot',plot.title = element_text(face='bold'))

liv.deRes.overlap.oppDir <- liv.deRes.overlap %>% 
  filter((log2FoldChange_PostDexExp * log2FoldChange_Jansen2019) < 0) # filter to genes with opposite directions 
# write_tsv(liv.deRes.overlap.oppDir,'analysis/gene_level_rnaseq/reversed_genes/liver_reversedGenes_02.28.23.tsv')

liv.deRes.overlap.oppDir%>% filter(log2FoldChange_PostDexExp > 0) %>% nrow()

liv.deRes.pdexpUnique <- liv.deRes.fed_vs_notFed %>% 
  filter(!(gene_number %in% liv.deRes.hib_vs_act$gene_number)) 
liv.deRes.jansenUnique <- liv.deRes.hib_vs_act %>% 
  filter(!(gene_number %in% liv.deRes.fed_vs_notFed$gene_number)) 


# Muscle ------------------------------------------------------------------

mus.deRes.fed_vs_notFed <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/muscle_PWRes_Fed.vs.NotFed.tsv')
mus.deRes.fed_vs_pdhib <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/muscle_PWRes_Fed.vs.PDHib.tsv')
mus.deRes.hib_vs_act <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/muscle_J19Hib.vs.J19Act.tsv')

mus.venn <- list('Post-Dex Fed vs.\nPost-Dex Not Fed' = mus.deRes.fed_vs_notFed$id,
                 # 'Post-Dex Fed vs.\nPost-Dex Hib' = mus.deRes.fed_vs_pdhib$id,
                 'Jansen 2019 Hib vs.\nJansen 2019 Act' = mus.deRes.hib_vs_act$id)

p.mus2 <- ggvenn(mus.venn,set_name_size = 4,fill_color=c('white','white','white'),show_percentage = F)+ ggtitle('Muscle') + theme(plot.title = element_text(face='bold',size=18))

mus.deRes.overlap <- mus.deRes.fed_vs_notFed %>% 
  filter(gene_number %in% mus.deRes.hib_vs_act$gene_number) %>% 
  left_join(mus.deRes.hib_vs_act,by=c('gene_number'),suffix = c('_PostDexExp','_Jansen2019'))

p.mus <- ggplot(mus.deRes.overlap,aes(y=log2FoldChange_PostDexExp,x=log2FoldChange_Jansen2019)) +
  geom_point(size=2,alpha=0.6,color='seagreen4') +
  geom_hline(yintercept = 0,lwd=0.2) + geom_vline(xintercept = 0,lwd=0.2)+
  geom_abline(slope = -1,lty=2,alpha=0.5)+  
  labs(title='',x='Active vs. Hibernation Log2FC',y=' ') +
  # coord_cartesian(xlim = c(-5,5),ylim=c(-3,3)) +
  theme_linedraw(base_size = 14) + theme(plot.title.position = 'plot',plot.title = element_text(face='bold'))
p.mus
p.adi2 + p.liv2 + p.mus2 + p.adi + p.liv + p.mus

mus.deRes.overlap.oppDir <- mus.deRes.overlap %>% 
  filter((log2FoldChange_PostDexExp * log2FoldChange_Jansen2019) < 0) # filter to genes with opposite directions 
# write_tsv(mus.deRes.overlap.oppDir,'analysis/gene_level_rnaseq/reversed_genes/muscle_reversedGenes_02.28.23.tsv')

mus.deRes.overlap.oppDir%>% filter(log2FoldChange_PostDexExp > 0) %>% nrow()

mus.deRes.pdexpUnique <- mus.deRes.fed_vs_notFed %>% 
  filter(!(gene_number %in% mus.deRes.hib_vs_act$gene_number)) 
mus.deRes.jansenUnique <- mus.deRes.hib_vs_act %>% 
  filter(!(gene_number %in% mus.deRes.fed_vs_notFed$gene_number)) 



# GO Analysis of Relevant Subsets -----------------------------------------

# List of Entrez IDs for conversion below
hs <- org.Hs.eg.db


# Define function for running GO, KEGG, Disease enrich --------------------

runAllFunctionalEnrich <- function(target_list,background_list) {
  # Run Biological process GO enrichment
  goRes.bp <- enrichGO(gene=target_list,
                       universe = background_list,
                       OrgDb = org.Hs.eg.db,
                       ont = 'BP',
                       pAdjustMethod = 'BH',
                       pvalueCutoff = 0.01,
                       qvalueCutoff = 0.05,
                       readable = T) %>% 
    simplify(cutoff=0.6, by="p.adjust", select_fun=min) %>% 
    mutate(category = 'Biological Process')
  
  # Run Molecular Function GO enrichment
  goRes.mf <- enrichGO(gene=target_list,
                       universe = background_list,
                       OrgDb = org.Hs.eg.db,
                       ont = 'MF',
                       pAdjustMethod = 'BH',
                       pvalueCutoff = 0.01,
                       qvalueCutoff = 0.05,
                       readable = T) %>% 
    simplify(cutoff=0.6, by="p.adjust", select_fun=min)%>% 
    mutate(category = 'Molecular Function')
  
  # Run Cellular Component GO enrichment
  goRes.cc <- enrichGO(gene=target_list,
                       universe = background_list,
                       OrgDb = org.Hs.eg.db,
                       ont = 'CC',
                       pAdjustMethod = 'BH',
                       pvalueCutoff = 0.01,
                       qvalueCutoff = 0.05,
                       readable = T) %>% 
    simplify(cutoff=0.6, by="p.adjust", select_fun=min) %>% 
    mutate(category = 'Cellular Component')
  
  # Run KEGG Enrichment
  kRes <- enrichKEGG(gene         = target_list,
                     universe = background_list,
                     organism     = 'hsa',
                     pvalueCutoff = 0.05)%>% 
    mutate(category = 'KEGG')
  
  # Run DisGeNet enrichment analysis
  dgnRes <- enrichDGN(target_list,universe = background_list) 
  
  # Combine GO results
  goRes.all <- data.frame(goRes.bp) %>% bind_rows(data.frame(goRes.mf),data.frame(goRes.cc),data.frame(kRes))
  
  out <- list('BP'=goRes.bp,'MF'=goRes.mf,'CC'=goRes.cc,'KEGG'=kRes,'DisGeNet' = dgnRes,'AllGO'=goRes.all)
  
  return(out)  
}




# Adipose GO analyses -----------------------------------------------------

adip.gene_all <- read_tsv('analysis/gene_level_rnaseq/norm_counts/adipose_vstNormCounts_08.11.22.tsv') %>% 
  mutate(gene_id = str_split_fixed(gene_id,'[:]',2)[,2])

# Make background gene list for GO analyses

adip.test.bg <- adip.gene_all %>% 
  select(gene_id) %>% 
  mutate(gene = str_split_fixed(gene_id,'[:]',2)[,2])

# Define genes of interest and background, convert to Entrez IDs
adip.bg <- adip.test.bg
adip.up.de <- adi.deRes.overlap.oppDir %>% filter(log2FoldChange_PostDexExp > 0)
adip.dw.de <- adi.deRes.overlap.oppDir %>% filter(log2FoldChange_PostDexExp < 0)


adip.bg.symbols <- adip.bg$gene

adip.bg.entrez <- AnnotationDbi::select(hs, 
                    keys = adip.bg.symbols,
                    columns = c("ENTREZID", "SYMBOL"),
                    keytype = "SYMBOL")


adip.up.de.symbols <- adip.up.de$symbol_PostDexExp
adip.up.de.entrez <- AnnotationDbi::select(hs, 
                     keys = adip.up.de.symbols,
                     columns = c("ENTREZID", "SYMBOL"),
                     keytype = "SYMBOL")

adip.dw.de.symbols <- adip.dw.de$symbol_PostDexExp
adip.dw.de.entrez <- AnnotationDbi::select(hs, 
                                           keys = adip.dw.de.symbols,
                                           columns = c("ENTREZID", "SYMBOL"),
                                           keytype = "SYMBOL")


# liver GO analyses -----------------------------------------------------

liv.gene_all <- read_tsv('analysis/gene_level_rnaseq/norm_counts/liver_vstNormCounts_08.11.22.tsv') %>% 
  mutate(gene_id = str_split_fixed(gene_id,'[:]',2)[,2])

# Make background gene list for GO analyses

liv.test.bg <- liv.gene_all %>% 
  select(gene_id) %>% 
  mutate(gene = str_split_fixed(gene_id,'[:]',2)[,2])

# Define genes of interest and background, convert to Entrez IDs
liv.bg <- liv.test.bg
liv.up.de <- liv.deRes.overlap.oppDir %>% filter(log2FoldChange_PostDexExp > 0)
liv.dw.de <- liv.deRes.overlap.oppDir %>% filter(log2FoldChange_PostDexExp < 0)

# write_tsv(liv.up.de,'analysis/gene_level_rnaseq/intersect_sets/liv_reversed_up.tsv')
# write_tsv(liv.dw.de,'analysis/gene_level_rnaseq/intersect_sets/liv_reversed_dw.tsv')


liv.bg.symbols <- liv.bg$gene
liv.bg.entrez <- AnnotationDbi::select(hs, 
                                        keys = liv.bg.symbols,
                                        columns = c("ENTREZID", "SYMBOL"),
                                        keytype = "SYMBOL")


liv.up.de.symbols <- liv.up.de$symbol_PostDexExp
liv.up.de.entrez <- AnnotationDbi::select(hs, 
                                           keys = liv.up.de.symbols,
                                           columns = c("ENTREZID", "SYMBOL"),
                                           keytype = "SYMBOL")

liv.dw.de.symbols <- liv.dw.de$symbol_PostDexExp
liv.dw.de.entrez <- AnnotationDbi::select(hs, 
                                           keys = liv.dw.de.symbols,
                                           columns = c("ENTREZID", "SYMBOL"),
                                           keytype = "SYMBOL")




# muscle GO analyses -----------------------------------------------------

mus.gene_all <- read_tsv('analysis/gene_level_rnaseq/norm_counts/muscle_vstNormCounts_08.11.22.tsv') %>% 
  mutate(gene_id = str_split_fixed(gene_id,'[:]',2)[,2])

# Make background gene list for GO analyses

mus.test.bg <- mus.gene_all %>% 
  select(gene_id) %>% 
  mutate(gene = str_split_fixed(gene_id,'[:]',2)[,2])

# Define genes of interest and background, convert to Entrez IDs
mus.bg <- mus.test.bg
mus.up.de <- mus.deRes.overlap.oppDir %>% filter(log2FoldChange_PostDexExp > 0)
mus.dw.de <- mus.deRes.overlap.oppDir %>% filter(log2FoldChange_PostDexExp < 0)

# write_tsv(mus.up.de,'analysis/gene_level_rnaseq/intersect_sets/mus_reversed_up.tsv')
# write_tsv(mus.dw.de,'analysis/gene_level_rnaseq/intersect_sets/mus_reversed_dw.tsv')

mus.bg.symbols <- mus.bg$gene
mus.bg.entrez <- AnnotationDbi::select(hs, 
                                       keys = mus.bg.symbols,
                                       columns = c("ENTREZID", "SYMBOL"),
                                       keytype = "SYMBOL")


mus.up.de.symbols <- mus.up.de$symbol_PostDexExp
mus.up.de.entrez <- AnnotationDbi::select(hs, 
                                          keys = mus.up.de.symbols,
                                          columns = c("ENTREZID", "SYMBOL"),
                                          keytype = "SYMBOL")

mus.dw.de.symbols <- mus.dw.de$symbol_PostDexExp
mus.dw.de.entrez <- AnnotationDbi::select(hs, 
                                          keys = mus.dw.de.symbols,
                                          columns = c("ENTREZID", "SYMBOL"),
                                          keytype = "SYMBOL")



### GO analysis with jansen 2019 genes as background

#Adipose
adip.allJansen.bg.symbol <- adi.deRes.hib_vs_act$symbol
adip.allJansen.bg.entrez <- AnnotationDbi::select(hs,
                                                 keys = adip.allJansen.bg.symbol,
                                                 columns = c("ENTREZID", "SYMBOL"),
                                                 keytype = "SYMBOL")


adip.up.JansenBG.funcRes <- runAllFunctionalEnrich(target_list = adip.up.de.entrez$ENTREZID,
                                                   background_list = adip.allJansen.bg.entrez$ENTREZID)
adip.dw.JansenBG.funcRes <- runAllFunctionalEnrich(target_list = adip.dw.de.entrez$ENTREZID,
                                                   background_list = adip.allJansen.bg.entrez$ENTREZID)


adi.combined <- adip.up.JansenBG.funcRes$AllGO %>%
  bind_rows(adip.dw.JansenBG.funcRes$AllGO,.id='id') %>%
  mutate(id = case_when(
    id == 1 ~ 'Upregulated Post-Feeding',
    id == 2 ~ 'Downregulated Post-Feeding'
  )) %>% 
  mutate(GeneRatio = sapply(GeneRatio, function(x) eval(parse(text=x))))%>% 
  mutate(category = factor(category,levels=c('Molecular Function','Biological Process','Cellular Component','KEGG')))


p.adi.combined <- ggplot(adi.combined, 
                         aes(x=GeneRatio, fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(fill=category),pch=21,size=3,show.legend = F) +
  scale_x_continuous(expand = c(0.02,0.02)) +
  labs(x="Gene Ratio",y=NULL) +
  ggtitle('A) Adipose') +
  facet_grid(cols = vars(id),rows=vars(category),scales = 'free_y',space = 'free_y') +
  theme_linedraw(base_size = 12) + theme(plot.title.position = 'plot')
p.adi.combined


#Liver
liv.allJansen.bg.symbol <- liv.deRes.hib_vs_act$symbol
liv.allJansen.bg.entrez <- AnnotationDbi::select(hs,
                                       keys = liv.allJansen.bg.symbol,
                                       columns = c("ENTREZID", "SYMBOL"),
                                       keytype = "SYMBOL")


liv.up.JansenBG.funcRes <- runAllFunctionalEnrich(target_list = liv.up.de.entrez$ENTREZID,
                                                  background_list = liv.allJansen.bg.entrez$ENTREZID)
liv.dw.JansenBG.funcRes <- runAllFunctionalEnrich(target_list = liv.dw.de.entrez$ENTREZID,
                                                  background_list = liv.allJansen.bg.entrez$ENTREZID)


liv.combined <- liv.up.JansenBG.funcRes$AllGO %>%
  bind_rows(liv.dw.JansenBG.funcRes$AllGO,.id='id') %>%
  mutate(id = case_when(
    id == 1 ~ 'Upregulated Post-Feeding',
    id == 2 ~ 'Downregulated Post-Feeding'
  )) %>% 
  mutate(GeneRatio = sapply(GeneRatio, function(x) eval(parse(text=x))))%>% 
  mutate(category = factor(category,levels=c('Molecular Function','Biological Process','Cellular Component','KEGG')))


p.liv.combined <- ggplot(liv.combined, 
                         aes(x=GeneRatio, fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(fill=category),pch=21,size=3,show.legend = F) +
  scale_x_continuous(expand = c(0.02,0.02)) +
  labs(x="Gene Ratio",y=NULL) +
  ggtitle('B) Liver') +
  facet_grid(cols = vars(id),rows=vars(category),scales = 'free_y',space = 'free_y') +
  theme_linedraw(base_size = 12) + theme(plot.title.position = 'plot')
p.liv.combined


#Muscle
mus.allJansen.bg.symbol <- mus.deRes.hib_vs_act$symbol
mus.allJansen.bg.entrez <- AnnotationDbi::select(hs,
                                                 keys = mus.allJansen.bg.symbol,
                                                 columns = c("ENTREZID", "SYMBOL"),
                                                 keytype = "SYMBOL")


mus.up.JansenBG.funcRes <- runAllFunctionalEnrich(target_list = mus.up.de.entrez$ENTREZID,
                                                  background_list = mus.allJansen.bg.entrez$ENTREZID)
mus.dw.JansenBG.funcRes <- runAllFunctionalEnrich(target_list = mus.dw.de.entrez$ENTREZID,
                                                  background_list = mus.allJansen.bg.entrez$ENTREZID)

mus.combined <- mus.up.JansenBG.funcRes$AllGO %>%
  bind_rows(mus.dw.JansenBG.funcRes$AllGO,.id='id') %>%
  mutate(id = case_when(
    id == 1 ~ 'Upregulated Post-Feeding',
    id == 2 ~ 'Downregulated Post-Feeding'
  )) %>% 
  mutate(GeneRatio = sapply(GeneRatio, function(x) eval(parse(text=x))))%>% 
  mutate(category = factor(category,levels=c('Molecular Function','Biological Process','Cellular Component','KEGG')))


p.mus.combined <- ggplot(mus.combined, 
                         aes(x=GeneRatio, fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(fill=category),pch=21,size=3,show.legend = F) +
  scale_x_continuous(expand = c(0.02,0.02)) +
  labs(x="Gene Ratio",y=NULL) +
  ggtitle('C) Muscle') +
  facet_grid(cols = vars(id),rows=vars(category),scales = 'free_y',space = 'free_y') +
  theme_linedraw(base_size = 12) + theme(plot.title.position = 'plot')
p.mus.combined


p1 <- p.adi.combined / p.mus.combined + plot_layout(heights = c(1,8))

p2 <- p.liv.combined 

p1 - p2 


#

