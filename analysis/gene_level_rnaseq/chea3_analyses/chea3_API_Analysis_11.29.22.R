library(httr)
library(jsonlite)
library(tidyverse)

# Adipose -----------------------------------------------------------------

adi.deRes.fed_vs_notFed <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/adipose_PWRes_Fed.vs.NotFed.tsv')
adi.deRes.hib_vs_act <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/adipose_J19Hib.vs.J19Act.tsv')

adi.deRes.overlap <- adi.deRes.fed_vs_notFed %>% 
  filter(gene_number %in% adi.deRes.hib_vs_act$gene_number) %>% 
  left_join(adi.deRes.hib_vs_act,by=c('gene_number'),suffix = c('_PostDexExp','_Jansen2019'))

adi.deRes.overlap.oppDir <- adi.deRes.overlap %>% 
  filter((log2FoldChange_PostDexExp * log2FoldChange_Jansen2019) < 0) # filter to genes with opposite directions 

adi.postfed.reverse.up <- adi.deRes.overlap.oppDir %>% filter(log2FoldChange_PostDexExp > 0) %>% select(symbol_PostDexExp)
adi.postfed.reverse.dw <- adi.deRes.overlap.oppDir %>% filter(log2FoldChange_PostDexExp < 0) %>% select(symbol_PostDexExp)

adi.hib.up <- adi.deRes.hib_vs_act %>% filter(log2FoldChange > 0) %>% select(symbol)
adi.hib.dw <- adi.deRes.hib_vs_act %>% filter(log2FoldChange < 0) %>% select(symbol)

# write_tsv(adi.postfed.reverse.up,'analysis/gene_level_rnaseq/chea3_analyses/input_lists/adi.reversed.up.txt',col_names = F)
# write_tsv(adi.postfed.reverse.dw,'analysis/gene_level_rnaseq/chea3_analyses/input_lists/adi.reversed.dw.txt',col_names = F)
# write_tsv(adi.hib.up,'analysis/gene_level_rnaseq/chea3_analyses/input_lists/adi.hib.up.txt',col_names = F)
# write_tsv(adi.hib.dw,'analysis/gene_level_rnaseq/chea3_analyses/input_lists/adi.hib.dw.txt',col_names = F)



# Liver -----------------------------------------------------------------

liv.deRes.fed_vs_notFed <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/liver_PWRes_Fed.vs.NotFed.tsv')
liv.deRes.hib_vs_act <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/liver_J19Hib.vs.J19Act.tsv')

liv.deRes.overlap <- liv.deRes.fed_vs_notFed %>% 
  filter(gene_number %in% liv.deRes.hib_vs_act$gene_number) %>% 
  left_join(liv.deRes.hib_vs_act,by=c('gene_number'),suffix = c('_PostDexExp','_Jansen2019'))

liv.deRes.overlap.oppDir <- liv.deRes.overlap %>% 
  filter((log2FoldChange_PostDexExp * log2FoldChange_Jansen2019) < 0) # filter to genes with opposite directions 

liv.postfed.reverse.up <- liv.deRes.overlap.oppDir %>% filter(log2FoldChange_PostDexExp > 0) %>% select(symbol_PostDexExp)
liv.postfed.reverse.dw <- liv.deRes.overlap.oppDir %>% filter(log2FoldChange_PostDexExp < 0) %>% select(symbol_PostDexExp)

liv.hib.up <- liv.deRes.hib_vs_act %>% filter(log2FoldChange > 0) %>% select(symbol)
liv.hib.dw <- liv.deRes.hib_vs_act %>% filter(log2FoldChange < 0) %>% select(symbol)

# write_tsv(liv.postfed.reverse.up,'analysis/gene_level_rnaseq/chea3_analyses/input_lists/liv.reversed.up.txt',col_names = F)
# write_tsv(liv.postfed.reverse.dw,'analysis/gene_level_rnaseq/chea3_analyses/input_lists/liv.reversed.dw.txt',col_names = F)
# write_tsv(liv.hib.up,'analysis/gene_level_rnaseq/chea3_analyses/input_lists/liv.hib.up.txt',col_names = F)
# write_tsv(liv.hib.dw,'analysis/gene_level_rnaseq/chea3_analyses/input_lists/liv.hib.dw.txt',col_names = F)




# Muscle -----------------------------------------------------------------

mus.deRes.fed_vs_notFed <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/muscle_PWRes_Fed.vs.NotFed.tsv')
mus.deRes.hib_vs_act <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/muscle_J19Hib.vs.J19Act.tsv')

mus.deRes.overlap <- mus.deRes.fed_vs_notFed %>% 
  filter(gene_number %in% mus.deRes.hib_vs_act$gene_number) %>% 
  left_join(mus.deRes.hib_vs_act,by=c('gene_number'),suffix = c('_PostDexExp','_Jansen2019'))

mus.deRes.overlap.oppDir <- mus.deRes.overlap %>% 
  filter((log2FoldChange_PostDexExp * log2FoldChange_Jansen2019) < 0) # filter to genes with opposite directions 

mus.postfed.reverse.up <- mus.deRes.overlap.oppDir %>% filter(log2FoldChange_PostDexExp > 0) %>% select(symbol_PostDexExp)
mus.postfed.reverse.dw <- mus.deRes.overlap.oppDir %>% filter(log2FoldChange_PostDexExp < 0) %>% select(symbol_PostDexExp)

mus.hib.up <- mus.deRes.hib_vs_act %>% filter(log2FoldChange > 0) %>% select(symbol)
mus.hib.dw <- mus.deRes.hib_vs_act %>% filter(log2FoldChange < 0) %>% select(symbol)

# write_tsv(mus.postfed.reverse.up,'analysis/gene_level_rnaseq/chea3_analyses/input_lists/mus.reversed.up.txt',col_names = F)
# write_tsv(mus.postfed.reverse.dw,'analysis/gene_level_rnaseq/chea3_analyses/input_lists/mus.reversed.dw.txt',col_names = F)
# write_tsv(mus.hib.up,'analysis/gene_level_rnaseq/chea3_analyses/input_lists/mus.hib.up.txt',col_names = F)
# write_tsv(mus.hib.dw,'analysis/gene_level_rnaseq/chea3_analyses/input_lists/mus.hib.dw.txt',col_names = F)




# CHEA3 API Analyses ------------------------------------------------------

# Define function to run CHEA3 using API and return mean_rank result table

runChea3 <- function(gene_list,query_name) {
  url = "https://maayanlab.cloud/chea3/api/enrich/"
  encode = "json"
  payload = list(query_name = query_name, gene_set = gene_list)
  response = POST(url = url, body = payload, encode = encode)
  json = content(response, "text")
  
  #results as list of R dataframes
  results = fromJSON(json)
  
  #mean rank results
  mean_rank = results$`Integrated--meanRank`
  return(mean_rank)
}

## Adipose
adi.reverse.up.chea3res <- runChea3(adi.postfed.reverse.up$symbol_PostDexExp,'adipose - post-fed reversed, up')
adi.reverse.dw.chea3res <- runChea3(adi.postfed.reverse.dw$symbol_PostDexExp,'adipose - post-fed reversed, down')
adi.hib.up.chea3res <- runChea3(adi.hib.up$symbol,'adipose - hibernation, up')
adi.hib.dw.chea3res <- runChea3(adi.hib.dw$symbol,'adipose - hibernation, down')

# write_tsv(adi.reverse.up.chea3res,'./analysis/gene_level_rnaseq/chea3_analyses/results/adi.reverse.up.chea3res.txt')
# write_tsv(adi.reverse.dw.chea3res,'./analysis/gene_level_rnaseq/chea3_analyses/results/adi.reverse.dw.chea3res.txt')
# write_tsv(adi.hib.up.chea3res,'./analysis/gene_level_rnaseq/chea3_analyses/results/adi.hib.up.chea3res.txt')
# write_tsv(adi.hib.dw.chea3res,'./analysis/gene_level_rnaseq/chea3_analyses/results/adi.hib.dw.chea3res.txt')


## Liver
liv.reverse.up.chea3res <- runChea3(liv.postfed.reverse.up$symbol_PostDexExp,'liver - post-fed reversed, up')
liv.reverse.dw.chea3res <- runChea3(liv.postfed.reverse.dw$symbol_PostDexExp,'liver - post-fed reversed, down')
liv.hib.up.chea3res <- runChea3(liv.hib.up$symbol,'liver - hibernation, up')
liv.hib.dw.chea3res <- runChea3(liv.hib.dw$symbol,'liver - hibernation, down')

# write_tsv(liv.reverse.up.chea3res,'./analysis/gene_level_rnaseq/chea3_analyses/results/liv.reverse.up.chea3res.txt')
# write_tsv(liv.reverse.dw.chea3res,'./analysis/gene_level_rnaseq/chea3_analyses/results/liv.reverse.dw.chea3res.txt')
# write_tsv(liv.hib.up.chea3res,'./analysis/gene_level_rnaseq/chea3_analyses/results/liv.hib.up.chea3res.txt')
# write_tsv(liv.hib.dw.chea3res,'./analysis/gene_level_rnaseq/chea3_analyses/results/liv.hib.dw.chea3res.txt')


## Muscle
mus.reverse.up.chea3res <- runChea3(mus.postfed.reverse.up$symbol_PostDexExp,'muscle - post-fed reversed, up')
mus.reverse.dw.chea3res <- runChea3(mus.postfed.reverse.dw$symbol_PostDexExp,'muscle - post-fed reversed, down')
mus.hib.up.chea3res <- runChea3(mus.hib.up$symbol,'muscle - hibernation, up')
mus.hib.dw.chea3res <- runChea3(mus.hib.dw$symbol,'muscle - hibernation, down')

# write_tsv(mus.reverse.up.chea3res,'./analysis/gene_level_rnaseq/chea3_analyses/results/mus.reverse.up.chea3res.txt')
# write_tsv(mus.reverse.dw.chea3res,'./analysis/gene_level_rnaseq/chea3_analyses/results/mus.reverse.dw.chea3res.txt')
# write_tsv(mus.hib.up.chea3res,'./analysis/gene_level_rnaseq/chea3_analyses/results/mus.hib.up.chea3res.txt')
# write_tsv(mus.hib.dw.chea3res,'./analysis/gene_level_rnaseq/chea3_analyses/results/mus.hib.dw.chea3res.txt')

## Candidate proteins

candprot.list <- c('IGF1', 'IGFALS', 'C1S', 'C2', 'SOD3', 'VTN', 'IGFBP2', 'JCHAIN')

candProt.chea3res <- runChea3(candprot.list,'candidate proteins')
# write_tsv(candProt.chea3res,'./analysis/gene_level_rnaseq/chea3_analyses/results/candidateProteins.chea3res.txt')
