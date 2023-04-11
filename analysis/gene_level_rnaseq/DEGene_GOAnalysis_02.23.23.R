library(tidyverse)
# install.packages("BiocManager")
# BiocManager::install("DOSE")
# BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
# devtools::install_github("YuLab-SMU/clusterProfiler")
library('clusterProfiler')
library(DOSE)
library(patchwork)


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

adip.de.fed_vs_notFed <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/adipose_PWRes_Fed.vs.NotFed.tsv')

adip.up.de <- adip.de.fed_vs_notFed %>% filter(log2FoldChange > 0)
adip.dw.de <- adip.de.fed_vs_notFed %>% filter(log2FoldChange < 0)


adip.bg.symbols <- adip.bg$gene
adip.bg.entrez <- AnnotationDbi::select(hs, 
                                        keys = adip.bg.symbols,
                                        columns = c("ENTREZID", "SYMBOL"),
                                        keytype = "SYMBOL")


adip.up.de.symbols <- adip.up.de$symbol
adip.up.de.entrez <- AnnotationDbi::select(hs, 
                                           keys = adip.up.de.symbols,
                                           columns = c("ENTREZID", "SYMBOL"),
                                           keytype = "SYMBOL")

adip.dw.de.symbols <- adip.dw.de$symbol
adip.dw.de.entrez <- AnnotationDbi::select(hs, 
                                           keys = adip.dw.de.symbols,
                                           columns = c("ENTREZID", "SYMBOL"),
                                           keytype = "SYMBOL")



# Run GO/functional enrich analyses ---------------------------------------------------------


adip.up.funcRes <- runAllFunctionalEnrich(adip.up.de.entrez$ENTREZID,adip.bg.entrez$ENTREZID)
adip.down.funcRes <- runAllFunctionalEnrich(adip.dw.de.entrez$ENTREZID,adip.bg.entrez$ENTREZID)


adip.combined <- adip.up.funcRes$AllGO %>% 
  bind_rows(adip.down.funcRes$AllGO,.id='id') %>% 
  mutate(id = case_when(
    id == 1 ~ 'Upregulated Post-Feeding',
    id == 2 ~ 'Downregulated Post-Feeding'
  )) %>% 
  mutate(GeneRatio = sapply(GeneRatio, function(x) eval(parse(text=x)))) %>% 
  mutate(category = factor(category,levels=c('Molecular Function','Biological Process','Cellular Component','KEGG')))

p.adi.combined <- ggplot(adip.combined, showCategory = 10, 
       aes(x=GeneRatio, fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(fill=category),pch=21,size=3,show.legend = F) +
  scale_x_continuous(expand = c(0.02,0.02)) +
  labs(x="Gene Ratio",y=NULL) +
  ggtitle('A) Adipose') +
  facet_grid(cols = vars(id),rows=vars(category),scales = 'free_y',space = 'free_y') +
  theme_linedraw(base_size = 12) + theme(plot.title.position = 'plot')

p.adi.combined


# liver GO analyses -----------------------------------------------------

liv.gene_all <- read_tsv('analysis/gene_level_rnaseq/norm_counts/liver_vstNormCounts_08.11.22.tsv') %>% 
  mutate(gene_id = str_split_fixed(gene_id,'[:]',2)[,2])

# Make background gene list for GO analyses

liv.test.bg <- liv.gene_all %>% 
  select(gene_id) %>% 
  mutate(gene = str_split_fixed(gene_id,'[:]',2)[,2])

# Define genes of interest and background, convert to Entrez IDs
liv.bg <- liv.test.bg

liv.de.fed_vs_notFed <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/liver_PWRes_Fed.vs.NotFed.tsv')

liv.up.de <- liv.de.fed_vs_notFed %>% filter(log2FoldChange > 0)
liv.dw.de <- liv.de.fed_vs_notFed %>% filter(log2FoldChange < 0)


liv.bg.symbols <- liv.bg$gene
liv.bg.entrez <- AnnotationDbi::select(hs, 
                                        keys = liv.bg.symbols,
                                        columns = c("ENTREZID", "SYMBOL"),
                                        keytype = "SYMBOL")


liv.up.de.symbols <- liv.up.de$symbol
liv.up.de.entrez <- AnnotationDbi::select(hs, 
                                           keys = liv.up.de.symbols,
                                           columns = c("ENTREZID", "SYMBOL"),
                                           keytype = "SYMBOL")

liv.dw.de.symbols <- liv.dw.de$symbol
liv.dw.de.entrez <- AnnotationDbi::select(hs, 
                                           keys = liv.dw.de.symbols,
                                           columns = c("ENTREZID", "SYMBOL"),
                                           keytype = "SYMBOL")



# Run GO/functional enrich analyses ---------------------------------------------------------

liv.up.funcRes <- runAllFunctionalEnrich(liv.up.de.entrez$ENTREZID,liv.bg.entrez$ENTREZID)
liv.down.funcRes <- runAllFunctionalEnrich(liv.dw.de.entrez$ENTREZID,liv.bg.entrez$ENTREZID)


liv.combined <- liv.up.funcRes$AllGO %>% 
  bind_rows(liv.down.funcRes$AllGO,.id='id') %>% 
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



# muscle GO analyses -----------------------------------------------------

mus.gene_all <- read_tsv('analysis/gene_level_rnaseq/norm_counts/muscle_vstNormCounts_08.11.22.tsv') %>% 
  mutate(gene_id = str_split_fixed(gene_id,'[:]',2)[,2])

# Make background gene list for GO analyses

mus.test.bg <- mus.gene_all %>% 
  select(gene_id) %>% 
  mutate(gene = str_split_fixed(gene_id,'[:]',2)[,2])

# Define genes of interest and background, convert to Entrez IDs
mus.bg <- mus.test.bg

mus.de.fed_vs_notFed <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/muscle_PWRes_Fed.vs.NotFed.tsv')

mus.up.de <- mus.de.fed_vs_notFed %>% filter(log2FoldChange > 0)
mus.dw.de <- mus.de.fed_vs_notFed %>% filter(log2FoldChange < 0)


mus.bg.symbols <- mus.bg$gene
mus.bg.entrez <- AnnotationDbi::select(hs, 
                                       keys = mus.bg.symbols,
                                       columns = c("ENTREZID", "SYMBOL"),
                                       keytype = "SYMBOL")


mus.up.de.symbols <- mus.up.de$symbol
mus.up.de.entrez <- AnnotationDbi::select(hs, 
                                          keys = mus.up.de.symbols,
                                          columns = c("ENTREZID", "SYMBOL"),
                                          keytype = "SYMBOL")

mus.dw.de.symbols <- mus.dw.de$symbol
mus.dw.de.entrez <- AnnotationDbi::select(hs, 
                                          keys = mus.dw.de.symbols,
                                          columns = c("ENTREZID", "SYMBOL"),
                                          keytype = "SYMBOL")



# Run GO/functional enrich analyses ---------------------------------------------------------

mus.up.funcRes <- runAllFunctionalEnrich(mus.up.de.entrez$ENTREZID,mus.bg.entrez$ENTREZID)
mus.down.funcRes <- runAllFunctionalEnrich(mus.dw.de.entrez$ENTREZID,mus.bg.entrez$ENTREZID)


mus.combined <- mus.up.funcRes$AllGO %>% 
  bind_rows(mus.down.funcRes$AllGO,.id='id') %>% 
  mutate(id = case_when(
    id == 1 ~ 'Upregulated Post-Feeding',
    id == 2 ~ 'Downregulated Post-Feeding'
  )) %>% 
  mutate(GeneRatio = sapply(GeneRatio, function(x) eval(parse(text=x))))%>% 
  mutate(category = factor(category,levels=c('Molecular Function','Biological Process','Cellular Component','KEGG')))


p.mus.combined <- ggplot(mus.combined, showCategory = 10, 
       aes(x=GeneRatio, fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(fill=category),pch=21,size=3,show.legend = F) +
  scale_x_continuous(expand = c(0.02,0.02)) +
  labs(x="Gene Ratio",y=NULL) +
  ggtitle('C) Muscle') +
  facet_grid(cols = vars(id),rows=vars(category),scales = 'free_y',space = 'free_y') +
  theme_linedraw(base_size = 12) + theme(plot.title.position = 'plot')



# Supp Fig plotting -------------------------------------------------------

p1 <- p.adi.combined / p.liv.combined + plot_layout(heights = c(1,8))

p2 <- p.mus.combined 

p1 - p2 





