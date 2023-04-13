library(tidyverse)
library(org.Hs.eg.db)
# BiocManager::install("clusterProfiler")
library('clusterProfiler')
library(DOSE)


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
    mutate(category = 'biological_process')
  
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
    mutate(category = 'molecular_function')
  
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
    mutate(category = 'cellular_component')
  
  # Run KEGG Enrichment
  kRes <- enrichKEGG(gene         = target_list,
                     universe = background_list,
                     organism     = 'hsa',
                     pvalueCutoff = 0.05)
  
  # Run DisGeNet enrichment analysis
  dgnRes <- enrichDGN(target_list,universe = background_list) 
  
  # Combine GO results
  goRes.all <- data.frame(goRes.bp) %>% bind_rows(data.frame(goRes.mf),data.frame(goRes.cc))
  
  out <- list('BP'=goRes.bp,'MF'=goRes.mf,'CC'=goRes.cc,'KEGG'=kRes,'DisGeNet' = dgnRes,'AllGO'=goRes.all)
  
  return(out)  
}


# Adipose GO analyses -----------------------------------------------------

adip.gene_all <- read_tsv('analysis/isoform_level_rnaseq/adipose.filteredIsoforms.tsv') %>% 
  mutate(gene = str_split_fixed(gene_id,'[-]',2)[,2])

# Make background gene list for GO analyses

adip.test.bg <- adip.gene_all %>% 
  select(gene) %>% 
  unique()

# Define genes of interest and background, convert to Entrez IDs
adip.bg <- adip.test.bg

adip.de.fed_vs_notFed <- read_tsv('analysis/isoform_level_rnaseq/isoformSwitchResults/adipose_allSigSwitches_02.24.23.tsv') %>% 
  filter(condition_1 == 'Fat_PDExp_Fed' & condition_2 == 'Fat_PDExp_NotFed') %>% 
  mutate(gene = str_split_fixed(gene_id,'[-]',2)[,2])


adip.bg.symbols <- adip.bg$gene
adip.bg.entrez <- AnnotationDbi::select(hs, 
                                        keys = adip.bg.symbols,
                                        columns = c("ENTREZID", "SYMBOL"),
                                        keytype = "SYMBOL")


adip.sig.de.symbols <- adip.de.fed_vs_notFed$gene
adip.sig.de.entrez <- AnnotationDbi::select(hs, 
                                           keys = adip.sig.de.symbols,
                                           columns = c("ENTREZID", "SYMBOL"),
                                           keytype = "SYMBOL")





# Run GO/functional enrich analyses ---------------------------------------------------------

adip.sig.funcRes <- runAllFunctionalEnrich(adip.sig.de.entrez$ENTREZID,adip.bg.entrez$ENTREZID)

adip.sig.funcRes

ggplot(adip.sig.funcRes$BP, showCategory = 15, 
       aes(x=GeneRatio, fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(fill=p.adjust, size = Count),shape=21) +
  scale_fill_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(5, 10)) +
  scale_x_continuous(expand = c(0.02,0.02)) +
  labs(x="Gene Ratio",y=NULL,title='Biological Process') +
  theme_linedraw(base_size = 16) + theme(plot.title.position = 'plot')

ggplot(adip.sig.funcRes$CC, showCategory = 15, 
       aes(x=GeneRatio, fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(fill=p.adjust, size = Count),shape=21) +
  scale_fill_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(5, 10)) +
  scale_x_continuous(expand = c(0.02,0.02)) +
  labs(x="Gene Ratio",y=NULL,title='Cellular Component') +
  theme_linedraw(base_size = 16) + theme(plot.title.position = 'plot')

ggplot(adip.sig.funcRes$MF, showCategory = 15, 
       aes(x=GeneRatio, fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(fill=p.adjust, size = Count),shape=21) +
  scale_fill_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(5, 10)) +
  scale_x_continuous(expand = c(0.02,0.02)) +
  labs(x="Gene Ratio",y=NULL,title='Molecular Function') +
  theme_linedraw(base_size = 16) + theme(plot.title.position = 'plot')


adip.sig.funcRes$KEGG

adip.sig.funcRes$DisGeNet


# liver GO analyses -----------------------------------------------------


liv.gene_all <- read_tsv('analysis/isoform_level_rnaseq/liver.filteredIsoforms.tsv') %>% 
  mutate(gene = str_split_fixed(gene_id,'[-]',2)[,2])

# Make background gene list for GO analyses

liv.test.bg <- liv.gene_all %>% 
  select(gene) %>% 
  unique()

# Define genes of interest and background, convert to Entrez IDs
liv.bg <- liv.test.bg

liv.de.fed_vs_notFed <- read_tsv('analysis/isoform_level_rnaseq/isoformSwitchResults/liver_allSigSwitches_02.24.23.tsv') %>% 
  filter(condition_1 == 'Liver_PDExp_Fed' & condition_2 == 'Liver_PDExp_NotFed') %>% 
  mutate(gene = str_split_fixed(gene_id,'[-]',2)[,2])


liv.bg.symbols <- liv.bg$gene
liv.bg.entrez <- AnnotationDbi::select(hs, 
                                        keys = liv.bg.symbols,
                                        columns = c("ENTREZID", "SYMBOL"),
                                        keytype = "SYMBOL")


liv.sig.de.symbols <- liv.de.fed_vs_notFed$gene
liv.sig.de.entrez <- AnnotationDbi::select(hs, 
                                            keys = liv.sig.de.symbols,
                                            columns = c("ENTREZID", "SYMBOL"),
                                            keytype = "SYMBOL")





# Run GO/functional enrich analyses ---------------------------------------------------------

liv.sig.funcRes <- runAllFunctionalEnrich(liv.sig.de.entrez$ENTREZID,liv.bg.entrez$ENTREZID)

liv.sig.funcRes

ggplot(liv.sig.funcRes$BP, showCategory = 15, 
       aes(x=GeneRatio, fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(fill=p.adjust, size = Count),shape=21) +
  scale_fill_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(5, 10)) +
  scale_x_continuous(expand = c(0.02,0.02)) +
  labs(x="Gene Ratio",y=NULL,title='Biological Process') +
  theme_linedraw(base_size = 16) + theme(plot.title.position = 'plot')

ggplot(liv.sig.funcRes$CC, showCategory = 15, 
       aes(x=GeneRatio, fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(fill=p.adjust, size = Count),shape=21) +
  scale_fill_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(5, 10)) +
  scale_x_continuous(expand = c(0.02,0.02)) +
  labs(x="Gene Ratio",y=NULL,title='Cellular Component') +
  theme_linedraw(base_size = 16) + theme(plot.title.position = 'plot')

ggplot(liv.sig.funcRes$MF, showCategory = 15, 
       aes(x=GeneRatio, fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(fill=p.adjust, size = Count),shape=21) +
  scale_fill_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(5, 10)) +
  scale_x_continuous(expand = c(0.02,0.02)) +
  labs(x="Gene Ratio",y=NULL,title='Molecular Function') +
  theme_linedraw(base_size = 16) + theme(plot.title.position = 'plot')


liv.sig.funcRes$KEGG

liv.sig.funcRes$DisGeNet




# muscle GO analyses -----------------------------------------------------



mus.gene_all <- read_tsv('analysis/isoform_level_rnaseq/muscle.filteredIsoforms.tsv') %>% 
  mutate(gene = str_split_fixed(gene_id,'[-]',2)[,2])

# Make background gene list for GO analyses

mus.test.bg <- mus.gene_all %>% 
  select(gene) %>% 
  unique()

# Define genes of interest and background, convert to Entrez IDs
mus.bg <- mus.test.bg

mus.de.fed_vs_notFed <- read_tsv('analysis/isoform_level_rnaseq/isoformSwitchResults/muscle_allSigSwitches_02.24.23.tsv') %>% 
  filter(condition_1 == 'Muscle_PDExp_Fed' & condition_2 == 'Muscle_PDExp_NotFed') %>% 
  mutate(gene = str_split_fixed(gene_id,'[-]',2)[,2])


mus.bg.symbols <- mus.bg$gene
mus.bg.entrez <- AnnotationDbi::select(hs, 
                                       keys = mus.bg.symbols,
                                       columns = c("ENTREZID", "SYMBOL"),
                                       keytype = "SYMBOL")


mus.sig.de.symbols <- mus.de.fed_vs_notFed$gene
mus.sig.de.entrez <- AnnotationDbi::select(hs, 
                                           keys = mus.sig.de.symbols,
                                           columns = c("ENTREZID", "SYMBOL"),
                                           keytype = "SYMBOL")





# Run GO/functional enrich analyses ---------------------------------------------------------

mus.sig.funcRes <- runAllFunctionalEnrich(mus.sig.de.entrez$ENTREZID,mus.bg.entrez$ENTREZID)

mus.sig.funcRes

ggplot(mus.sig.funcRes$BP, showCategory = 15, 
       aes(x=GeneRatio, fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(fill=p.adjust, size = Count),shape=21) +
  scale_fill_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(5, 10)) +
  scale_x_continuous(expand = c(0.02,0.02)) +
  labs(x="Gene Ratio",y=NULL,title='Biological Process') +
  theme_linedraw(base_size = 16) + theme(plot.title.position = 'plot')

ggplot(mus.sig.funcRes$CC, showCategory = 15, 
       aes(x=GeneRatio, fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(fill=p.adjust, size = Count),shape=21) +
  scale_fill_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(5, 10)) +
  scale_x_continuous(expand = c(0.02,0.02)) +
  labs(x="Gene Ratio",y=NULL,title='Cellular Component') +
  theme_linedraw(base_size = 16) + theme(plot.title.position = 'plot')

ggplot(mus.sig.funcRes$MF, showCategory = 15, 
       aes(x=GeneRatio, fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(fill=p.adjust, size = Count),shape=21) +
  scale_fill_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(5, 10)) +
  scale_x_continuous(expand = c(0.02,0.02)) +
  labs(x="Gene Ratio",y=NULL,title='Molecular Function') +
  theme_linedraw(base_size = 16) + theme(plot.title.position = 'plot')


mus.sig.funcRes$KEGG

mus.sig.funcRes$DisGeNet


