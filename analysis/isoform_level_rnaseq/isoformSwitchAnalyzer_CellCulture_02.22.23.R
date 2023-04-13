
library(tidyverse)
library(IsoformSwitchAnalyzeR)
library(tximport)
library(gt)

# Read in Kallisto quantification files ------------------------------------
txi.kallisto <- importIsoformExpression(parentDir = 'data/4_kallisto/cellCulture_kallisto_quants/',addIsofomIdAsColumn = T)

txi.kallisto$counts <- txi.kallisto$counts
txi.kallisto$abundance <- txi.kallisto$abundance 
txi.kallisto$length <- txi.kallisto$length 

sample_info <- txi.kallisto$counts %>% 
  pivot_longer(-1,names_to = 'sample',values_to = 'count') %>% 
  dplyr::select(2) %>% 
  unique() %>% 
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
  ))


# Build design matrix -------------------------------------------

designMat <- sample_info %>% 
  dplyr::select(sample,treatment,sex) %>% 
  dplyr::select(sampleID=1,condition=2,sex)

comparisons <- t(combn(unique(designMat$condition),m=2)) %>% 
  as.data.frame() %>% 
  dplyr::select(condition_1=1,condition_2=2)



# Build switch list object ------------------------------------------------
aSwitchList <- importRdata(
  isoformCountMatrix   = txi.kallisto$counts,
  isoformRepExpression = txi.kallisto$abundance,
  designMatrix         = designMat,
  isoformExonAnnoation = 'data/0_reference/GCF_023065955.1_UrsArc1.0_genomic.gtf',
  isoformNtFasta       = 'data/4_kallisto/GCF_023065955.1_UrsArc1.0_rna.fa',
  showProgress = T,
  addAnnotatedORFs = T,
  comparisonsToMake = comparisons
)

# Pre-filtering -----------------------------------------------------------
aSwitchListFiltered <- preFilter(
  switchAnalyzeRlist = aSwitchList,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE
)


# Perform isoform switch tests using DEXseq method ------------------------
aSwitchListAnalyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = aSwitchListFiltered,
  reduceToSwitchingGenes=TRUE
)
switchSummary <- extractSwitchSummary(aSwitchListAnalyzed)




# Extract amino acid sequences per isoform --------------------------------

aSwitchListAnalyzed <- extractSequence(
  aSwitchListAnalyzed,
  pathToOutput = 'analysis/isoform_level_rnaseq/isoform_seqs',
  writeToFile=T,
  # writeToFile=F # to avoid output when running this example data

)


# Incorporate PFAM results ------------------------------------------------

# Ran pfam with this command:
# pfam_scan.pl -fasta ../isoform_seqs/isoformSwitchAnalyzeR_isoform_AA.fasta -dir . -outfile ./cellCulture_pfam_scan_res_02.22.23.txt

aSwitchListAnalyzed <- analyzePFAM(
  switchAnalyzeRlist   = aSwitchListAnalyzed,
  pathToPFAMresultFile = 'analysis/isoform_level_rnaseq/cellCulture_pfam_scan_res_02.22.23.txt',
  showProgress=T
)


# Predicting Alternative Splicing -----------------------------------------

aSwitchListAnalyzed <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = aSwitchListAnalyzed,
  quiet=F
)


## overview of number of intron retentions (IR)
table(aSwitchListAnalyzed$AlternativeSplicingAnalysis$IR )


# Predicting Switch Consequences ------------------------------------------

consequencesOfInterest <- c('intron_retention','NMD_status','domains_identified','ORF_seq_similarity',
                            'tss','tts',
                            '5_utr_length',
                            '3_utr_length'
                            )

aSwitchListAnalyzed <- analyzeSwitchConsequences(
  aSwitchListAnalyzed,
  consequencesToAnalyze = consequencesOfInterest,
  dIFcutoff = 0.1,
  showProgress=T
)


extractSwitchSummary(aSwitchListAnalyzed, dIFcutoff = 0.1, filterForConsequences = FALSE)
extractSwitchSummary(aSwitchListAnalyzed, dIFcutoff = 0.1, filterForConsequences = TRUE)


# Adipose: Analysis of Individual Isoform Switching --------------------------------

# Subset to only adipose Post-Dex treatments
adipose.SwitchListAnalyzedSubset <- subsetSwitchAnalyzeRlist(
  aSwitchListAnalyzed,
  str_detect(aSwitchListAnalyzed$isoformFeatures$condition_1,'HH')
)


# Extract all genes with switches, sorted by Qval
adipose.PDExp.allSigSwitches <- extractTopSwitches(
  adipose.SwitchListAnalyzedSubset,
  filterForConsequences = TRUE,
  n = NA,
  sortByQvals = T
)
write_tsv(adipose.PDExp.allSigSwitches,'analysis/isoform_level_rnaseq/isoformSwitchResults/cellCulture_allSigSwitches_02.22.23.tsv')


save(designMat, comparisons, aSwitchList, aSwitchListFiltered, aSwitchListAnalyzed, adipose.SwitchListAnalyzedSubset,adipose.PDExp.allSigSwitches, switchSummary, file = "analysis/isoform_level_rnaseq/cellCulture_IsoformSwitchAnalyzer.RData")
# load('analysis/isoform_level_rnaseq/cellCulture_IsoformSwitchAnalyzer.RData')

adipose.PDExp.allSigSwitches.summary <- adipose.PDExp.allSigSwitches %>% 
  # filter(str_detect(condition_1,'50Per',negate = T) & str_detect(condition_2,'50Per',negate = T)) %>% 
  group_by(condition_1,condition_2) %>% 
  tally()  %>% 
  dplyr::select(1,2,`of consequence`=n)

adipose.PDExp.allSigSwitches.notConsFiltered <- extractTopSwitches(
  adipose.SwitchListAnalyzedSubset,
  filterForConsequences = F,
  n = NA,
  sortByQvals = T
)

adipose.PDExp.allSigSwitches.notConsFiltered.summary <- adipose.PDExp.allSigSwitches.notConsFiltered %>% 
  filter(str_detect(condition_1,'50Per',negate = T) & str_detect(condition_2,'50Per',negate = T)) %>% 
  group_by(condition_1,condition_2) %>% 
  tally() %>% 
  dplyr::select(1,2,all=n)

adipose.PDExp.allSigSwitches.both.summary <- adipose.PDExp.allSigSwitches.summary %>% 
  left_join(adipose.PDExp.allSigSwitches.notConsFiltered.summary,by=c('condition_1','condition_2'))

adipose.PDExp.allSigSwitches.both.summary %>% 
  ungroup() %>%  
  dplyr::select(1,2,4,3) %>% 
  mutate(condition_1 = str_remove_all(condition_1,'Fat_'),
         condition_2 = str_remove_all(condition_2,'Fat_')) %>% 
  gt() %>% 
  tab_header(title = 'Adipose') %>% 
  tab_spanner(label = "# of Genes with Isoform Switches", columns = matches("all|of consequence")) %>% 
  tab_source_note(md("*Of consequence*: change in intron retention, NMD status, protein domains, tss, tts, and/or UTR length"))



# Comparison to tissue-level

adipose.PDExp.allSigSwitches

tissue.PDExp.allSigSwitches <- read_tsv('analysis/isoform_level_rnaseq/isoformSwitchResults/adipose_allSigSwitches_02.24.23.tsv') %>% 
  filter(condition_2 == 'Fat_PDExp_NotFed' & condition_1 == 'Fat_PDExp_Fed')


overlap.switches <- adipose.PDExp.allSigSwitches %>% 
  filter(gene_id %in% tissue.PDExp.allSigSwitches$gene_id)

