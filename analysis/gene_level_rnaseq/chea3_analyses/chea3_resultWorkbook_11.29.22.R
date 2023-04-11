
library(tidyverse)
library(pheatmap)
library(network)
library(ggnetwork)
library(colorspace)
library(ggvenn)
library(patchwork)

# Define top N number of TFs to consider ----------------------------------
top_n = 25

# Adipose -----------------------------------------------------------------


adi.hib.up.tfs <- read_tsv('analysis/gene_level_rnaseq/chea3_analyses/results/adi.hib.up.chea3res.txt') %>% janitor::clean_names() %>% filter(rank <= top_n)
adi.hib.dw.tfs <- read_tsv('analysis/gene_level_rnaseq/chea3_analyses/results/adi.hib.dw.chea3res.txt') %>% janitor::clean_names() %>% filter(rank <= top_n)

adi.fed.up.tfs <- read_tsv('analysis/gene_level_rnaseq/chea3_analyses/results/adi.reverse.up.chea3res.txt') %>% janitor::clean_names() %>% filter(rank <= top_n)
adi.fed.dw.tfs <- read_tsv('analysis/gene_level_rnaseq/chea3_analyses/results/adi.reverse.dw.chea3res.txt') %>% janitor::clean_names() %>% filter(rank <= top_n)

adi.fed.both.tfs <- adi.fed.up.tfs %>% bind_rows(adi.fed.dw.tfs)


# Chea3 tissue-associated TFs
tissueEnrichDB <- read_tsv('analysis/gene_level_rnaseq/chea3_analyses/databases/wgcna_gtex_annotated5.reformat.txt')
tissueEnrichDB.adi <- tissueEnrichDB %>% filter(tissue_general=='Adipose Tissue')

# sum(adi.hib.up.tfs$tf %in% tissueEnrichDB.adi$tf) # 0 / 25
# sum(adi.hib.dw.tfs$tf %in% tissueEnrichDB.adi$tf) # 2 / 25
sum(adi.fed.up.tfs$tf %in% tissueEnrichDB.adi$tf) # 7 / 25
sum(adi.fed.dw.tfs$tf %in% tissueEnrichDB.adi$tf) # 0 / 25


# Liver -----------------------------------------------------------------

liv.hib.up.tfs <- read_tsv('analysis/gene_level_rnaseq/chea3_analyses/results/liv.hib.up.chea3res.txt') %>% janitor::clean_names() %>% filter(rank <= top_n)
liv.hib.dw.tfs <- read_tsv('analysis/gene_level_rnaseq/chea3_analyses/results/liv.hib.dw.chea3res.txt') %>% janitor::clean_names() %>% filter(rank <= top_n)

liv.fed.up.tfs <- read_tsv('analysis/gene_level_rnaseq/chea3_analyses/results/liv.reverse.up.chea3res.txt') %>% janitor::clean_names() %>% filter(rank <= top_n)
liv.fed.dw.tfs <- read_tsv('analysis/gene_level_rnaseq/chea3_analyses/results/liv.reverse.dw.chea3res.txt') %>% janitor::clean_names() %>% filter(rank <= top_n)

liv.fed.both.tfs <- liv.fed.up.tfs %>% bind_rows(liv.fed.dw.tfs)


# Chea3 tissue-associated TFs
tissueEnrichDB.liv <- tissueEnrichDB %>% filter(tissue_general=='Liver')

# sum(liv.hib.up.tfs$tf %in% tissueEnrichDB.liv$tf) # 0 / 25
# sum(liv.hib.dw.tfs$tf %in% tissueEnrichDB.liv$tf) # 10 / 25
sum(liv.fed.up.tfs$tf %in% tissueEnrichDB.liv$tf) # 10 / 25
sum(liv.fed.dw.tfs$tf %in% tissueEnrichDB.liv$tf) # 4 / 25


# Muscle -----------------------------------------------------------------

mus.hib.up.tfs <- read_tsv('analysis/gene_level_rnaseq/chea3_analyses/results/mus.hib.up.chea3res.txt') %>% janitor::clean_names() %>% filter(rank <= top_n)
mus.hib.dw.tfs <- read_tsv('analysis/gene_level_rnaseq/chea3_analyses/results/mus.hib.dw.chea3res.txt') %>% janitor::clean_names() %>% filter(rank <= top_n)

mus.fed.up.tfs <- read_tsv('analysis/gene_level_rnaseq/chea3_analyses/results/mus.reverse.up.chea3res.txt') %>% janitor::clean_names() %>% filter(rank <= top_n)
mus.fed.dw.tfs <- read_tsv('analysis/gene_level_rnaseq/chea3_analyses/results/mus.reverse.dw.chea3res.txt') %>% janitor::clean_names() %>% filter(rank <= top_n)

mus.fed.both.tfs <- mus.fed.up.tfs %>% bind_rows(mus.fed.dw.tfs)

# Chea3 tissue-associated TFs
tissueEnrichDB.mus <- tissueEnrichDB %>% filter(tissue_general=='Muscle')

# sum(mus.hib.up.tfs$tf %in% tissueEnrichDB.mus$tf) # 13 / 25
# sum(mus.hib.dw.tfs$tf %in% tissueEnrichDB.mus$tf) # 11 / 25
sum(mus.fed.up.tfs$tf %in% tissueEnrichDB.mus$tf) # 11 / 25
sum(mus.fed.dw.tfs$tf %in% tissueEnrichDB.mus$tf) # 8 / 25



# Overlap between tissues -------------------------------------------------

tf.venn.both <- list('Adipose' = adi.fed.both.tfs$tf,
                 'Liver' = liv.fed.both.tfs$tf,
                'Muscle'= mus.fed.both.tfs$tf)

ggvenn(tf.venn.both,set_name_size = 4,fill_color=c('white','white','white'),show_percentage = F) + theme(plot.title = element_text(face='bold',size=18))


tf.venn.up <- list('Adipose' = adi.fed.up.tfs$tf,
                     'Liver' = liv.fed.up.tfs$tf,
                     'Muscle'= mus.fed.up.tfs$tf)

up.venn <- ggvenn(tf.venn.up,set_name_size = 4,fill_color=c('white','white','white'),show_percentage = F) + theme(plot.title = element_text(face='bold',size=18))

#alternative to venn
adi.fed.up.tfs %>% 
  bind_rows(liv.fed.up.tfs,mus.fed.up.tfs) %>% 
  mutate(tissue = str_split_fixed(query_name,' ',2)[,1] %>% str_to_sentence()) %>% 
  select(tf,tissue) %>% 
  ggplot(aes(x=tissue,y=tf)) +
  geom_point() +
  theme_linedraw()

adi.fed.up.tfs %>% filter(tf %in% liv.fed.up.tfs$tf & tf %in% mus.fed.up.tfs$tf) ##PPARG shared in all three
adi.fed.up.tfs %>% filter(tf %in% liv.fed.up.tfs$tf) ## AR and CEBPA in adipose and liver
adi.fed.up.tfs %>% filter(tf %in% mus.fed.up.tfs$tf) ##TBX18 and PRRX1 in adipose and muscle

tf.venn.dw <- list('Adipose' = adi.fed.dw.tfs$tf,
                   'Liver' = liv.fed.dw.tfs$tf,
                   'Muscle'= mus.fed.dw.tfs$tf)

ggvenn(tf.venn.dw,set_name_size = 4,fill_color=c('white','white','white'),show_percentage = F) + theme(plot.title = element_text(face='bold',size=18))

adi.fed.dw.tfs %>% filter(tf %in% liv.fed.dw.tfs$tf) ## AR and CEBPA in adipose and liver

### Supp Table 4

all.fed.both.tfs <- adi.fed.both.tfs %>% bind_rows(liv.fed.both.tfs,mus.fed.both.tfs) %>% 
  mutate(direction = str_split_fixed(query_name,' ',3)[,3]) %>% 
  mutate(tissue = str_split_fixed(query_name,' ',2)[,1]) %>% 
  select(tissue,direction,everything()) %>% 
  select(-query_name)

# write_tsv(all.fed.both.tfs,'analysis/gene_level_rnaseq/chea3_analyses/SuppTable4Data_AllReversedTFs.tsv')

# PPARG Exploration -------------------------------------------------------

adi.fed.both.tfs %>% bind_rows(liv.fed.both.tfs,mus.fed.both.tfs) %>% 
  filter(tf == 'PPARG') %>% 
  mutate(gene_count = str_count(overlapping_genes,',')+1)

pparg.info <- adi.fed.both.tfs %>% bind_rows(liv.fed.both.tfs,mus.fed.both.tfs) %>% 
  filter(tf == 'PPARG') %>% 
  separate_rows(overlapping_genes,sep='[,]') 

pparg.info %>% 
  group_by(overlapping_genes) %>% 
  tally() %>% 
  arrange(-n)


# Recreating TF co-regulation network plotting from CHEA3 online ----------
# Notes on how this works:
# - Edges that were supported by evidence from two or more different libraries were retained in the network. 
# - Edges are directed where ChIP-seq evidence supports the interaction and are undirected in the case of co-occurrence or co-expression evidence only.

# read_db <- function(path) {
#   df <- read_delim(path,delim = ' ',col_names = c('unsplit')) %>%
#     mutate(tf = str_split_fixed(unsplit,'\t',2)[,1],
#            targets = str_split_fixed(unsplit,'\t',2)[,2]) %>%
#     select(tf,targets) %>%
#     mutate(targets = str_replace_all(targets,'\t',',')) %>%
#     separate_rows(targets,sep=',') %>%
#     mutate(tf = ifelse(str_detect(tf,'_'),str_split_fixed(tf,'[_]',2)[,1],tf)) %>% 
#     unique()
#   return(df)
# }
# 
#### WILL NEED TO UN GZIP THESE
# archs4.db <- read_db('analysis/gene_level_rnaseq/chea3_analyses/databases/ARCHS4_Coexpression.gmt') %>% mutate(db = 'archs4', type = 'coexp')
# encode.db <- read_db('analysis/gene_level_rnaseq/chea3_analyses/databases/ENCODE_ChIP-seq.gmt') %>% mutate(db = 'encode', type = 'chip')
# enrichr.db <- read_db('analysis/gene_level_rnaseq/chea3_analyses/databases/Enrichr_Queries.gmt') %>% mutate(db = 'enrichr', type = 'other')
# gtex.db <- read_db('analysis/gene_level_rnaseq/chea3_analyses/databases/GTEx_Coexpression.gmt') %>% mutate(db = 'gtex', type = 'coexp')
# literature.db <- read_db('analysis/gene_level_rnaseq/chea3_analyses/databases/Literature_ChIP-seq.gmt') %>% mutate(db = 'lit', type = 'chip')
# remap.db <- read_db('analysis/gene_level_rnaseq/chea3_analyses/databases/ReMap_ChIP-seq.gmt') %>% mutate(db = 'remap', type = 'chip')
# 
# all.db <- archs4.db %>% bind_rows(encode.db,enrichr.db,gtex.db,literature.db,remap.db)
# write_tsv(all.db,'analysis/gene_level_rnaseq/chea3_analyses/databases/all_combined.chea3db.txt')

#### WILL NEED TO UN GZIP THIS FILE FIRST
all.db <- read_tsv('analysis/gene_level_rnaseq/chea3_analyses/databases/all_combined.chea3db.txt')

unique(all.db$db)


# Build network plotting function -----------------------------------------

plot_coregNet <- function(df) {
  net.df <- all.db %>% 
    filter(tf %in% df$tf & targets %in% df$tf) %>% # filter to interactions involving focal TFs
    mutate(tf_target_pair = ifelse(tf < targets, paste(tf,targets,sep='_'),paste(targets,tf,sep='_'))) %>%  # combine tf-target columns in alphabetic order
    select(tf,targets,tf_target_pair,db,type) %>% 
    group_by(tf_target_pair) %>% # group by tf-target pairs (direction doesn't matter here)
    mutate(count = n_distinct(db)) %>%  # add column with number of databases supporting each tf-target pair
    ungroup() %>% # ungroup by tf-target pairs column
    group_by(tf,targets) %>% # group by tf and target columns separately (direction matters now)
    dplyr::mutate(chip_yes = ifelse(type == 'chip',1,0)) %>%  # add column indicating whether at least one db has chip evidence for a directed interaction
    mutate(chip_count = sum(chip_yes)) %>% select(-chip_yes) %>% # get the sum of dbs with chip evidence for each directed interaction
    select(-db,-type) %>% # remove db and type columns
    unique() %>% # make sure all rows are unique
    ungroup() %>% # ungroup
    filter(tf != targets) %>% # filter out self interactions if present
    filter(count > 1) %>%  # filter to only interactions in at least two databases
    mutate(chip_arrow = ifelse(chip_count > 0, 1,0)) %>% # make column indicating whether there is at least one chip interaction (used to draw arrows)
    select(tf,targets,count,chip_arrow)
  
  net <- network(net.df)
  net.v2 <- ggnetwork(net)
  
  ggplot(net.v2,aes(x=x,y=y,xend=xend,yend=yend)) +
    geom_edges(data=net.v2[net.v2$chip_arrow == 0,],alpha=1,lwd=.7) +
    geom_edges(data=net.v2[net.v2$chip_arrow == 1,],color='dodgerblue',alpha=1,lwd=.7,arrow = arrow(length = unit(6, "pt"), type = "closed")) +
    # scale_color_continuous_sequential(palette = 'Viridis',rev=F,) +
    geom_nodes(size=4) +
    geom_nodes(data=net.v2 %>% filter(vertex.names=='PPARG'),color='firebrick',size=4.2) +
    geom_nodelabel_repel(aes( label = vertex.names),box.padding = unit(0.5, "lines"),size=3) +
    theme_blank()
}

adi.up.net <- plot_coregNet(adi.fed.up.tfs) + ggtitle('Adipose - Upregulated Post-feeding') 
adi.dw.net <- plot_coregNet(adi.fed.dw.tfs) + ggtitle('Adipose - Downregulated Post-feeding')

adi.up.net + adi.dw.net

liv.up.net <- plot_coregNet(liv.fed.up.tfs) + ggtitle('Liver - Upregulated Post-feeding')
liv.dw.net <- plot_coregNet(liv.fed.dw.tfs) + ggtitle('Liver - Downregulated Post-feeding')

liv.up.net + liv.dw.net

mus.up.net <- plot_coregNet(mus.fed.up.tfs) + ggtitle('Muscle - Upregulated Post-feeding')
mus.dw.net <- plot_coregNet(mus.fed.dw.tfs) + ggtitle('Muscle - Downregulated Post-feeding')

mus.up.net + mus.dw.net

(adi.up.net + adi.dw.net) / (liv.up.net + liv.dw.net) / (mus.up.net + mus.dw.net)

(up.venn + adi.up.net) / (liv.up.net + mus.up.net)

adi.dw.net + liv.dw.net + mus.dw.net


###########################################################################################################################################################################################


# Candidate protein CHEA3 results -----------------------------------------

db.key <- all.db %>% select(db,type) %>% unique()

candprot.list <- c('IGF1', 'IGFALS', 'C1S', 'C2', 'SOD3', 'VTN', 'IGFBP2', 'JCHAIN')

candprot.tfs <- read_tsv('analysis/gene_level_rnaseq/chea3_analyses/results/candidateProteins.chea3res.txt') %>% janitor::clean_names() %>% filter(rank <= top_n) %>% 
  left_join(tissueEnrichDB)
# plot_coregNet(candprot.tfs)

candprot.tfs.long <- candprot.tfs %>% 
  select(tf,library,overlapping_genes) %>% 
  separate_rows(overlapping_genes,sep=',') %>% 
  separate_rows(library,sep=';') %>% 
  mutate(db = case_when(
    str_detect(library,'ARCHS4') ~ 'archs4',
    str_detect(library,'ENCODE') ~ 'encode',
    str_detect(library,'Enrichr') ~ 'enrichr',
    str_detect(library,'GTEx') ~ 'gtex',
    str_detect(library,'Liter') ~ 'lit',
    str_detect(library,'ReMap') ~ 'remap',
  )) %>% 
  left_join(db.key) %>% 
  select(tf,targets=overlapping_genes,db,type) %>% 
  group_by(tf,targets) %>% 
  dplyr::mutate(count = n()) %>% 
  dplyr::mutate(chip_yes = ifelse(type == 'chip',1,0)) %>%  # add column indicating whether at least one db has chip evidence for a directed interaction
  mutate(chip_count = sum(chip_yes)) %>% select(-chip_yes) %>%  # get the sum of dbs with chip evidence for each directed interaction
  select(-db,-type) %>% # remove db and type columns
  unique() %>% # make sure all rows are unique
  ungroup() %>% # ungroup
  filter(tf != targets) %>% # filter out self interactions if present
  filter(count > 1) %>%  # filter to only interactions in at least two databases
  mutate(chip_arrow = ifelse(chip_count > 0, 1,0)) %>% # make column indicating whether there is at least one chip interaction (used to draw arrows)
  select(tf,targets,count,chip_arrow)

candprot.tfs.long


prot_net.df <- all.db %>%
  filter((tf %in% candprot.tfs$tf & targets %in% candprot.tfs$tf)) %>% # filter to interactions involving focal TFs 
  mutate(tf_target_pair = ifelse(tf < targets, paste(tf,targets,sep='_'),paste(targets,tf,sep='_'))) %>%  # combine tf-target columns in alphabetic order
  select(tf,targets,tf_target_pair,db,type) %>%
  group_by(tf_target_pair) %>% # group by tf-target pairs (direction doesn't matter here)
  mutate(count = n_distinct(db)) %>%  # add column with number of databases supporting each tf-target pair
  ungroup() %>% # ungroup by tf-target pairs column
  group_by(tf,targets) %>% # group by tf and target columns separately (direction matters now)
  dplyr::mutate(chip_yes = ifelse(type == 'chip',1,0)) %>%  # add column indicating whether at least one db has chip evidence for a directed interaction
  mutate(chip_count = sum(chip_yes)) %>% select(-chip_yes) %>% # get the sum of dbs with chip evidence for each directed interaction
  select(-db,-type) %>% # remove db and type columns
  unique() %>% # make sure all rows are unique
  ungroup() %>% # ungroup
  filter(tf != targets) %>% # filter out self interactions if present
  filter(count > 1) %>%  # filter to only interactions in at least two databases
  mutate(chip_arrow = ifelse(chip_count > 0, 1,0)) %>% # make column indicating whether there is at least one chip interaction (used to draw arrows)
  select(tf,targets,count,chip_arrow) %>% 
  bind_rows(candprot.tfs.long) %>%  # add TF -> candidate protein interactions
  left_join(tissueEnrichDB)

# Plotting tests

prot_net <- network(prot_net.df)
prot_net.v2 <- ggnetwork(prot_net,layout = "fruchtermanreingold",niter=10000) %>% 
  mutate(candidate = ifelse(vertex.names %in% candprot.list,'Candidate Protein','TF'))

ggplot(prot_net.v2,aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_edges(data=prot_net.v2[prot_net.v2$chip_arrow == 0,],alpha=1,lwd=0.5,color='grey50') +
  geom_edges(data=prot_net.v2[prot_net.v2$chip_arrow == 1,],alpha=1,lwd=1,arrow = arrow(length = unit(5, "pt"), type = "closed")) +
  # scale_color_continuous_sequential(palette = 'Viridis') +
  geom_nodes(aes(fill=candidate,shape=candidate),size=6) +
  scale_shape_manual(values=c('Candidate Protein'=22,'TF'=21)) +
  geom_nodelabel_repel(aes( label = vertex.names),box.padding = unit(0.5, "lines"),size=3) +
  theme_blank()


prot_net.v2 %>% select(vertex.names,candidate) %>%
  unique() %>%
  left_join(tissueEnrichDB,by=c('vertex.names'='tf')) %>%
  select(vertex.names,candidate,tissue_general) %>%
  write_tsv('analysis/gene_level_rnaseq/chea3_analyses/CandidateProtein_TFNetwork_NodeInfo.tsv')
write_tsv(prot_net.df,'analysis/gene_level_rnaseq/chea3_analyses/CandidateProtein_TFNetwork_EdgeInfo.tsv')


# Overlap with reversed gene TFs ------------------------------------------

protGene.overlap <- candprot.tfs %>% select(-GO_enrichment,-tissue_specific,-tissue_general) %>% 
  bind_rows(adi.fed.both.tfs,liv.fed.both.tfs,mus.fed.both.tfs) %>% 
  group_by(tf) %>% 
  add_count() %>% 
  arrange(-n)

protGene.overlap %>% filter(str_detect(query_name,'proteins')) %>% filter(n>1) 

protGene.overlap %>% 
  filter(n>1) %>% 
  filter(tf %in% candprot.tfs$tf) %>% 
  ggplot(aes(y=tf,x=query_name)) +
  geom_point(size=6) +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

tissue.simple <- adi.fed.both.tfs %>% 
  bind_rows(liv.fed.both.tfs,mus.fed.both.tfs) %>% 
  mutate(targets = paste('Reversed genes in ', str_split_fixed(query_name,' ',2)[,1],sep = '')) %>% 
  select(targets,tf) %>% 
  unique() %>% 
  filter(tf %in% candprot.tfs$tf) %>% 
  mutate(count = 1, chip_arrow = 3, GO_enrichment = 'NA', tissue_general = 'NA',tissue_specific = 'NA')


prot_net.df.wTissues <- prot_net.df %>% bind_rows(tissue.simple)
#  write_tsv(prot_net.df.wTissues,'analysis/gene_level_rnaseq/chea3_analyses/CandidateProtein_TFNetwork_EdgeInfo_wTissues.tsv')



###########################################################################################################################################################################################


# Gene Expression of TFs -----------------------------------------

# Adipose
adi.deRes.fed_vs_notFed.up <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/adipose_PWRes_Fed.vs.NotFed.tsv') %>% filter(log2FoldChange > 0)
adi.deRes.fed_vs_notFed.dw <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/adipose_PWRes_Fed.vs.NotFed.tsv') %>% filter(log2FoldChange < 0)

adi.deRes.hib_vs_act.up <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/adipose_J19Hib.vs.J19Act.tsv') %>% filter(log2FoldChange > 0)
adi.deRes.hib_vs_act.dw <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/adipose_J19Hib.vs.J19Act.tsv') %>% filter(log2FoldChange < 0)

sum(adi.hib.up.tfs$tf %in% adi.deRes.hib_vs_act.up$symbol)
sum(adi.hib.dw.tfs$tf %in% adi.deRes.hib_vs_act.dw$symbol)

sum(adi.fed.up.tfs$tf %in% adi.deRes.fed_vs_notFed.up$symbol)
sum(adi.fed.dw.tfs$tf %in% adi.deRes.fed_vs_notFed.dw$symbol)

# Liver
liv.deRes.fed_vs_notFed.up <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/liver_PWRes_Fed.vs.NotFed.tsv') %>% filter(log2FoldChange > 0)
liv.deRes.fed_vs_notFed.dw <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/liver_PWRes_Fed.vs.NotFed.tsv') %>% filter(log2FoldChange < 0)

liv.deRes.hib_vs_act.up <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/liver_J19Hib.vs.J19Act.tsv') %>% filter(log2FoldChange > 0)
liv.deRes.hib_vs_act.dw <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/liver_J19Hib.vs.J19Act.tsv') %>% filter(log2FoldChange < 0)

sum(liv.hib.up.tfs$tf %in% liv.deRes.hib_vs_act.up$symbol)
sum(liv.hib.dw.tfs$tf %in% liv.deRes.hib_vs_act.dw$symbol)

sum(liv.fed.up.tfs$tf %in% liv.deRes.fed_vs_notFed.up$symbol)
sum(liv.fed.dw.tfs$tf %in% liv.deRes.fed_vs_notFed.dw$symbol)

# Muscle
mus.deRes.fed_vs_notFed.up <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/muscle_PWRes_Fed.vs.NotFed.tsv') %>% filter(log2FoldChange > 0)
mus.deRes.fed_vs_notFed.dw <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/muscle_PWRes_Fed.vs.NotFed.tsv') %>% filter(log2FoldChange < 0)

mus.deRes.hib_vs_act.up <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/muscle_J19Hib.vs.J19Act.tsv') %>% filter(log2FoldChange > 0)
mus.deRes.hib_vs_act.dw <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/muscle_J19Hib.vs.J19Act.tsv') %>% filter(log2FoldChange < 0)

sum(mus.hib.up.tfs$tf %in% mus.deRes.hib_vs_act.up$symbol)
sum(mus.hib.dw.tfs$tf %in% mus.deRes.hib_vs_act.dw$symbol)

sum(mus.fed.up.tfs$tf %in% mus.deRes.fed_vs_notFed.up$symbol)
sum(mus.fed.dw.tfs$tf %in% mus.deRes.fed_vs_notFed.dw$symbol)

mus.deRes.fed_vs_notFed.up %>% filter(symbol %in% mus.fed.up.tfs$tf)
mus.deRes.fed_vs_notFed.dw %>% filter(symbol %in% mus.fed.dw.tfs$tf)


# Proteins
adi.deRes.fed_vs_notFed.up$symbol %in% candprot.tfs.long$tf %>% sum()
adi.deRes.fed_vs_notFed.dw$symbol %in% candprot.tfs.long$tf %>% sum()

liv.deRes.fed_vs_notFed.up$symbol %in% candprot.tfs.long$tf %>% sum()
liv.deRes.fed_vs_notFed.dw$symbol %in% candprot.tfs.long$tf %>% sum()

mus.deRes.fed_vs_notFed.up$symbol %in% candprot.tfs.long$tf %>% sum() # AEBP1
mus.deRes.fed_vs_notFed.up %>% filter(symbol %in% candprot.tfs.long$tf)


mus.deRes.fed_vs_notFed.dw$symbol %in% candprot.tfs.long$tf %>% sum()



# 



