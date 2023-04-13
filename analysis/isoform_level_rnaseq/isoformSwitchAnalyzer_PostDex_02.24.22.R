
library(tidyverse)
library(IsoformSwitchAnalyzeR)
library(tximport)
library(gt)

# Read in Kallisto quantification files ------------------------------------
txi.kallisto <- importIsoformExpression(parentDir = 'data/4_kallisto/kallisto_quants/',addIsofomIdAsColumn = T)

# Remove samples with fewer than 3M reads mapped, known contaminants, known outlier samples, and CFA (previously shown to contain hair follicle)
# and rename column names to more informative IDs

sample_info <- read_tsv('sampleInfo_08.11.22.txt')

filter_out <- c('CFA','D2_S26','H2_S57','E2_S34','D7_S31','H7_S62', # CFA (contaminated) and low mapping
                'C2_S18','E3_S35','F2_S42', # contaminated samples
                'F3_S43','H6_S61','A2_S2' # outliers in gene expression pcas
                )

txi.kallisto$counts <- txi.kallisto$counts %>% 
  pivot_longer(-1,names_to = 'orig_name',values_to = 'values') %>% 
  filter(!(orig_name %in% filter_out)) %>% 
  left_join(sample_info) %>% 
  dplyr::select(isoform_id,combined,values) %>% 
  pivot_wider(names_from = combined,values_from = values)

txi.kallisto$abundance <- txi.kallisto$abundance %>% 
  pivot_longer(-1,names_to = 'orig_name',values_to = 'values') %>% 
  filter(!(orig_name %in% filter_out)) %>% 
  left_join(sample_info) %>% 
  dplyr::select(isoform_id,combined,values) %>% 
  pivot_wider(names_from = combined,values_from = values)

txi.kallisto$length <- txi.kallisto$length %>% 
  pivot_longer(-1,names_to = 'orig_name',values_to = 'values') %>% 
  filter(!(orig_name %in% filter_out)) %>% 
  left_join(sample_info) %>% 
  dplyr::select(isoform_id,combined,values) %>% 
  pivot_wider(names_from = combined,values_from = values)


# Build separate design matrix for each tissue -------------------------------------------

adip.designMat <- sample_info %>%
  filter(tissue=='Fat') %>%
  filter(str_detect(treatment,'50Per',negate = T)) %>%
  mutate(combined_condition = paste(tissue,treatment,sep='_')) %>%
  dplyr::select(combined,combined_condition,sex) %>%
  dplyr::select(sampleID=1,condition=2,sex)

liv.designMat <- sample_info %>%
  filter(tissue=='Liver') %>%
  mutate(combined_condition = paste(tissue,treatment,sep='_')) %>%
  dplyr::select(combined,combined_condition,sex) %>%
  dplyr::select(sampleID=1,condition=2,sex)

mus.designMat <- sample_info %>%
  filter(tissue=='Muscle') %>%
  mutate(combined_condition = paste(tissue,treatment,sep='_')) %>%
  dplyr::select(combined,combined_condition,sex) %>%
  dplyr::select(sampleID=1,condition=2,sex)


adip.comparisons <- t(combn(unique(adip.designMat$condition),m=2)) %>%
  as.data.frame() %>%
  dplyr::select(condition_1=1,condition_2=2) %>%
  mutate(temp_1 = str_split_fixed(condition_1,'_',2)[,1],
         temp_2 = str_split_fixed(condition_2,'_',2)[,1]) %>%
  filter(temp_1 == temp_2) %>%
  dplyr::select(-temp_1,-temp_2) %>%
  filter((str_detect(condition_1,'Jansen') & str_detect(condition_2,'Jansen')) | (str_detect(condition_1,'PDExp') & str_detect(condition_2,'PDExp')))

liv.comparisons <- t(combn(unique(liv.designMat$condition),m=2)) %>%
  as.data.frame() %>%
  dplyr::select(condition_1=1,condition_2=2) %>%
  mutate(temp_1 = str_split_fixed(condition_1,'_',2)[,1],
         temp_2 = str_split_fixed(condition_2,'_',2)[,1]) %>%
  filter(temp_1 == temp_2) %>%
  dplyr::select(-temp_1,-temp_2) %>%
  filter((str_detect(condition_1,'Jansen') & str_detect(condition_2,'Jansen')) | (str_detect(condition_1,'PDExp') & str_detect(condition_2,'PDExp')))

mus.comparisons <- t(combn(unique(mus.designMat$condition),m=2)) %>%
  as.data.frame() %>%
  dplyr::select(condition_1=1,condition_2=2) %>%
  mutate(temp_1 = str_split_fixed(condition_1,'_',2)[,1],
         temp_2 = str_split_fixed(condition_2,'_',2)[,1]) %>%
  filter(temp_1 == temp_2) %>%
  dplyr::select(-temp_1,-temp_2) %>%
  filter((str_detect(condition_1,'Jansen') & str_detect(condition_2,'Jansen')) | (str_detect(condition_1,'PDExp') & str_detect(condition_2,'PDExp')))



# Build switch list objects ------------------------------------------------
adip.SwitchList <- importRdata(
  isoformCountMatrix   = txi.kallisto$counts,
  isoformRepExpression = txi.kallisto$abundance,
  designMatrix         = adip.designMat,
  isoformExonAnnoation = 'data/0_reference/GCF_023065955.1_UrsArc1.0_genomic.gtf',
  isoformNtFasta       = 'data/4_kallisto/GCF_023065955.1_UrsArc1.0_rna.fa',
  showProgress = T,
  addAnnotatedORFs = T,
  comparisonsToMake = adip.comparisons
)

liv.SwitchList <- importRdata(
  isoformCountMatrix   = txi.kallisto$counts,
  isoformRepExpression = txi.kallisto$abundance,
  designMatrix         = liv.designMat,
  isoformExonAnnoation = 'data/0_reference/GCF_023065955.1_UrsArc1.0_genomic.gtf',
  isoformNtFasta       = 'data/4_kallisto/GCF_023065955.1_UrsArc1.0_rna.fa',
  showProgress = T,
  addAnnotatedORFs = T,
  comparisonsToMake = liv.comparisons
)

mus.SwitchList <- importRdata(
  isoformCountMatrix   = txi.kallisto$counts,
  isoformRepExpression = txi.kallisto$abundance,
  designMatrix         = mus.designMat,
  isoformExonAnnoation = 'data/0_reference/GCF_023065955.1_UrsArc1.0_genomic.gtf',
  isoformNtFasta       = 'data/4_kallisto/GCF_023065955.1_UrsArc1.0_rna.fa',
  showProgress = T,
  addAnnotatedORFs = T,
  comparisonsToMake = mus.comparisons
)


# Pre-filtering -----------------------------------------------------------
adip.SwitchListFiltered <- preFilter(
  switchAnalyzeRlist = adip.SwitchList,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE
)

liv.SwitchListFiltered <- preFilter(
  switchAnalyzeRlist = liv.SwitchList,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE
)

mus.SwitchListFiltered <- preFilter(
  switchAnalyzeRlist = mus.SwitchList,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE
)


## Export filtered gene list to use as background in GO analyses

write_tsv(adip.SwitchListFiltered$isoformFeatures,'analysis/isoform_level_rnaseq/adipose.filteredIsoforms.tsv')
write_tsv(liv.SwitchListFiltered$isoformFeatures,'analysis/isoform_level_rnaseq/liver.filteredIsoforms.tsv')
write_tsv(mus.SwitchListFiltered$isoformFeatures,'analysis/isoform_level_rnaseq/muscle.filteredIsoforms.tsv')



# Perform isoform switch tests using DEXseq method ------------------------
adip.SwitchListAnalyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = adip.SwitchListFiltered,
  reduceToSwitchingGenes=TRUE
)
adip.switchSummary <- extractSwitchSummary(adip.SwitchListAnalyzed)

liv.SwitchListAnalyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = liv.SwitchListFiltered,
  reduceToSwitchingGenes=TRUE
)
liv.switchSummary <- extractSwitchSummary(liv.SwitchListAnalyzed)

mus.SwitchListAnalyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = mus.SwitchListFiltered,
  reduceToSwitchingGenes=TRUE
)
mus.switchSummary <- extractSwitchSummary(mus.SwitchListAnalyzed)


# Extract amino acid sequences per isoform --------------------------------

adip.SwitchListAnalyzed <- extractSequence(
  adip.SwitchListAnalyzed,
  pathToOutput = 'analysis/isoform_level_rnaseq/isoform_seqs',
  writeToFile=T,
  outputPrefix='adipose'
  # writeToFile=F # to avoid output when running this example data
)

liv.SwitchListAnalyzed <- extractSequence(
  liv.SwitchListAnalyzed,
  pathToOutput = 'analysis/isoform_level_rnaseq/isoform_seqs',
  writeToFile=T,
  outputPrefix='liver'
  # writeToFile=F # to avoid output when running this example data
)

mus.SwitchListAnalyzed <- extractSequence(
  mus.SwitchListAnalyzed,
  pathToOutput = 'analysis/isoform_level_rnaseq/isoform_seqs',
  writeToFile=T,
  outputPrefix='muscle'
  # writeToFile=F # to avoid output when running this example data
)


# Incorporate PFAM results ------------------------------------------------

# Ran pfam with this command:
# pfam_scan.pl -fasta ../isoform_seqs/adipose_AA.fasta -dir . -outfile ./adipose_pfam_scan_res_02.24.23.txt
# pfam_scan.pl -fasta ../isoform_seqs/liver_AA.fasta -dir . -outfile ./liver_pfam_scan_res_02.24.23.txt
# pfam_scan.pl -fasta ../isoform_seqs/muscle_AA.fasta -dir . -outfile ./muscle_pfam_scan_res_02.24.23.txt
#

adip.SwitchListAnalyzed <- analyzePFAM(
  switchAnalyzeRlist   = adip.SwitchListAnalyzed,
  pathToPFAMresultFile = 'analysis/isoform_level_rnaseq/pfam_database/adipose_pfam_scan_res_02.24.23.txt',
  showProgress=T
)

liv.SwitchListAnalyzed <- analyzePFAM(
  switchAnalyzeRlist   = liv.SwitchListAnalyzed,
  pathToPFAMresultFile = 'analysis/isoform_level_rnaseq/pfam_database/liver_pfam_scan_res_02.24.23.txt',
  showProgress=T
)

mus.SwitchListAnalyzed <- analyzePFAM(
  switchAnalyzeRlist   = mus.SwitchListAnalyzed,
  pathToPFAMresultFile = 'analysis/isoform_level_rnaseq/pfam_database/muscle_pfam_scan_res_02.24.23.txt',
  showProgress=T
)



# Predicting Alternative Splicing -----------------------------------------

adip.SwitchListAnalyzed <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = adip.SwitchListAnalyzed,
  quiet=F
)

liv.SwitchListAnalyzed <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = liv.SwitchListAnalyzed,
  quiet=F
)

mus.SwitchListAnalyzed <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = mus.SwitchListAnalyzed,
  quiet=F
)


## overview of number of intron retentions (IR)
table(adip.SwitchListAnalyzed$AlternativeSplicingAnalysis$IR )
table(liv.SwitchListAnalyzed$AlternativeSplicingAnalysis$IR )
table(mus.SwitchListAnalyzed$AlternativeSplicingAnalysis$IR )


# Predicting Switch Consequences ------------------------------------------

consequencesOfInterest <- c('intron_retention','NMD_status','domains_identified','ORF_seq_similarity',
                            'tss','tts',
                            '5_utr_length',
                            '3_utr_length'
                            )

adip.SwitchListAnalyzed <- analyzeSwitchConsequences(
  adip.SwitchListAnalyzed,
  consequencesToAnalyze = consequencesOfInterest,
  dIFcutoff = 0.1,
  showProgress=T
)

liv.SwitchListAnalyzed <- analyzeSwitchConsequences(
  liv.SwitchListAnalyzed,
  consequencesToAnalyze = consequencesOfInterest,
  dIFcutoff = 0.1,
  showProgress=T
)

mus.SwitchListAnalyzed <- analyzeSwitchConsequences(
  mus.SwitchListAnalyzed,
  consequencesToAnalyze = consequencesOfInterest,
  dIFcutoff = 0.1,
  showProgress=T
)


extractSwitchSummary(adip.SwitchListAnalyzed, dIFcutoff = 0.1, filterForConsequences = FALSE)
extractSwitchSummary(adip.SwitchListAnalyzed, dIFcutoff = 0.1, filterForConsequences = TRUE)

extractSwitchSummary(liv.SwitchListAnalyzed, dIFcutoff = 0.1, filterForConsequences = FALSE)
extractSwitchSummary(liv.SwitchListAnalyzed, dIFcutoff = 0.1, filterForConsequences = TRUE)

extractSwitchSummary(mus.SwitchListAnalyzed, dIFcutoff = 0.1, filterForConsequences = FALSE)
extractSwitchSummary(mus.SwitchListAnalyzed, dIFcutoff = 0.1, filterForConsequences = TRUE)

save(adip.designMat, liv.designMat, mus.designMat, adip.comparisons, liv.comparisons, mus.comparisons, adip.SwitchList, liv.SwitchList, mus.SwitchList, adip.SwitchListFiltered, liv.SwitchListFiltered, mus.SwitchListFiltered, adip.SwitchListAnalyzed, liv.SwitchListAnalyzed, mus.SwitchListAnalyzed, adip.switchSummary, liv.switchSummary, mus.switchSummary, file = "analysis/isoform_level_rnaseq/IsoformSwitchAnalyzer.RData")
# load('analysis/isoform_level_rnaseq/IsoformSwitchAnalyzer.RData')


# Adipose: Analysis of Individual Isoform Switching --------------------------------

# Subset to only adipose Post-Dex treatments
adipose.SwitchListAnalyzedSubset <- adip.SwitchListAnalyzed

# Extract all genes with switches, sorted by Qval
adipose.PDExp.allSigSwitches <- extractTopSwitches(
  adipose.SwitchListAnalyzedSubset,
  filterForConsequences = TRUE,
  n = NA,
  sortByQvals = T
)
# write_tsv(adipose.PDExp.allSigSwitches,'analysis/isoform_level_rnaseq/isoformSwitchResults/adipose_allSigSwitches_02.24.23.tsv')


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

# plot isoform for insulin related gene
switchPlot(adipose.SwitchListAnalyzedSubset, gene = 'gene-GSK3B',condition2 = 'Fat_PDExp_Fed',condition1 = 'Fat_PDExp_NotFed')

cons <- adipose.SwitchListAnalyzedSubset$switchConsequence




# Liver: Analysis of Individual Isoform Switching --------------------------------

# Subset to only liver Post-Dex treatments
liver.SwitchListAnalyzedSubset <- liv.SwitchListAnalyzed
  
# Extract all genes with switches, sorted by Qval
liver.PDExp.allSigSwitches <- extractTopSwitches(
  liver.SwitchListAnalyzedSubset,
  filterForConsequences = TRUE,
  n = NA,
  sortByQvals = T
)
# write_tsv(liver.PDExp.allSigSwitches,'analysis/isoform_level_rnaseq/isoformSwitchResults/liver_allSigSwitches_02.24.23.tsv')

liver.PDExp.allSigSwitches.summary <- liver.PDExp.allSigSwitches %>% 
  # filter(str_detect(condition_1,'50Per',negate = T) & str_detect(condition_2,'50Per',negate = T)) %>% 
  group_by(condition_1,condition_2) %>% 
  tally()  %>% 
  dplyr::select(1,2,`of consequence`=n)

liver.PDExp.allSigSwitches.notConsFiltered <- extractTopSwitches(
  liver.SwitchListAnalyzedSubset,
  filterForConsequences = F,
  n = NA,
  sortByQvals = T
)

liver.PDExp.allSigSwitches.notConsFiltered.summary <- liver.PDExp.allSigSwitches.notConsFiltered %>% 
  filter(str_detect(condition_1,'50Per',negate = T) & str_detect(condition_2,'50Per',negate = T)) %>% 
  group_by(condition_1,condition_2) %>% 
  tally() %>% 
  dplyr::select(1,2,all=n)

liver.PDExp.allSigSwitches.both.summary <- liver.PDExp.allSigSwitches.summary %>% 
  left_join(liver.PDExp.allSigSwitches.notConsFiltered.summary,by=c('condition_1','condition_2'))

liver.PDExp.allSigSwitches.both.summary %>% 
  ungroup() %>%  
  dplyr::select(1,2,4,3) %>% 
  mutate(condition_1 = str_remove_all(condition_1,'Liver_'),
         condition_2 = str_remove_all(condition_2,'Liver_')) %>% 
  gt() %>% 
  tab_header(title = 'Liver') %>% 
  tab_spanner(label = "# of Genes with Isoform Switches", columns = matches("all|of consequence")) %>% 
  tab_source_note(md("*Of consequence*: change in intron retention, NMD status, protein domains, tss, tts, and/or UTR length"))



# muscle: Analysis of Individual Isoform Switching --------------------------------

# Subset to only muscle Post-Dex treatments
muscle.SwitchListAnalyzedSubset <- mus.SwitchListAnalyzed

# Extract all genes with switches, sorted by Qval
muscle.PDExp.allSigSwitches <- extractTopSwitches(
  muscle.SwitchListAnalyzedSubset,
  filterForConsequences = TRUE,
  n = NA,
  sortByQvals = T
)

# write_tsv(muscle.PDExp.allSigSwitches,'analysis/isoform_level_rnaseq/isoformSwitchResults/muscle_allSigSwitches_02.24.23.tsv')


muscle.PDExp.allSigSwitches.summary <- muscle.PDExp.allSigSwitches %>% 
  # filter(str_detect(condition_1,'50Per',negate = T) & str_detect(condition_2,'50Per',negate = T)) %>% 
  group_by(condition_1,condition_2) %>% 
  tally()  %>% 
  dplyr::select(1,2,`of consequence`=n)

muscle.PDExp.allSigSwitches.notConsFiltered <- extractTopSwitches(
  muscle.SwitchListAnalyzedSubset,
  filterForConsequences = F,
  n = NA,
  sortByQvals = T
)

muscle.PDExp.allSigSwitches.notConsFiltered.summary <- muscle.PDExp.allSigSwitches.notConsFiltered %>% 
  filter(str_detect(condition_1,'50Per',negate = T) & str_detect(condition_2,'50Per',negate = T)) %>% 
  group_by(condition_1,condition_2) %>% 
  tally() %>% 
  dplyr::select(1,2,all=n)

muscle.PDExp.allSigSwitches.both.summary <- muscle.PDExp.allSigSwitches.summary %>% 
  left_join(muscle.PDExp.allSigSwitches.notConsFiltered.summary,by=c('condition_1','condition_2'))

muscle.PDExp.allSigSwitches.both.summary %>% 
  ungroup() %>%  
  dplyr::select(1,2,4,3) %>% 
  mutate(condition_1 = str_remove_all(condition_1,'Muscle_'),
         condition_2 = str_remove_all(condition_2,'Muscle_')) %>% 
  gt() %>% 
  tab_header(title = 'muscle') %>% 
  tab_spanner(label = "# of Genes with Isoform Switches", columns = matches("all|of consequence")) %>% 
  tab_source_note(md("*Of consequence*: change in intron retention, NMD status, protein domains, tss, tts, and/or UTR length"))


# Combined result table for all tissues -----------------------------------

allTissue.PDExp.allSigSwitches.summary <- 
  muscle.PDExp.allSigSwitches.both.summary %>% 
  bind_rows(adipose.PDExp.allSigSwitches.both.summary,liver.PDExp.allSigSwitches.both.summary) %>% 
  mutate(tissue = str_split_fixed(condition_1,'[_]',2)[,1]) %>% 
  mutate(condition_1 = str_split_fixed(condition_1,'[_]',3)[,3],
         condition_2 = str_split_fixed(condition_2,'[_]',3)[,3]) %>% 
  ungroup() %>% 
  dplyr::select(tissue,`Condition 1`=condition_1,`Condition 2`=condition_2,3,4)

summary_table <- allTissue.PDExp.allSigSwitches.summary %>% 
  group_by(tissue) %>% 
  gt() %>% 
  tab_spanner(label = "# of Genes with Isoform Switches", columns = matches("all|of consequence")) %>% 
  tab_source_note(md("*Of consequence*: change in intron retention, NMD status, protein domains, tss, tts, and/or UTR length"))

summary_table

# summary_table %>% gtsave("analysis/isoform_level_rnaseq/AllTissue_IsoformSwitchSummary.html")



# Genome wide summaries of switches ---------------------------------------
## ADIPOSE
extractConsequenceSummary(
  adipose.SwitchListAnalyzedSubset,
  asFractionTotal = FALSE,
  plotGenes=FALSE,
)
extractConsequenceEnrichment(adipose.SwitchListAnalyzedSubset) # Nothing significant in feeding experiment

extractSwitchOverlap(
  subsetSwitchAnalyzeRlist(adipose.SwitchListAnalyzedSubset,
                           str_detect(adipose.SwitchListAnalyzedSubset$isoformFeatures$condition_1,'Jansen2019_Hib|PDExp_Fed') & 
                          str_detect(adipose.SwitchListAnalyzedSubset$isoformFeatures$condition_2,'Jansen2019_Act|PDExp_NotFed')),
  filterForConsequences=TRUE,
  plotIsoforms = T
)

adipose.overlap <- adipose.PDExp.allSigSwitches %>% 
  filter(switchConsequencesGene==TRUE & gene_switch_q_value < 0.05) %>% 
  filter((condition_1 == 'Fat_Jansen2019_Hib' & condition_2 == 'Fat_Jansen2019_Act') | (condition_1=='Fat_PDExp_Fed' & condition_2 == 'Fat_PDExp_NotFed')) %>% 
  mutate(comparison = paste(condition_1,':',condition_2,sep = '')) %>% 
  group_by(gene_id) %>% 
  add_count() %>% 
  arrange(-n)

adipose.overlap %>% filter(n==2)


## Liver
extractConsequenceSummary(
  liver.SwitchListAnalyzedSubset,
  asFractionTotal = FALSE,
  plotGenes=FALSE,
)
extractConsequenceEnrichment(liver.SwitchListAnalyzedSubset) # Nothing significant 

extractSwitchOverlap(
  subsetSwitchAnalyzeRlist(liver.SwitchListAnalyzedSubset,
                           str_detect(liver.SwitchListAnalyzedSubset$isoformFeatures$condition_1,'Jansen2019_Act|PDExp_Fed') & 
                             str_detect(liver.SwitchListAnalyzedSubset$isoformFeatures$condition_2,'Jansen2019_Hib|PDExp_NotFed')),
  filterForConsequences=TRUE,
  plotIsoforms = T
)

liver.overlap <- liver.PDExp.allSigSwitches %>% 
  filter(switchConsequencesGene==TRUE & gene_switch_q_value < 0.05) %>% 
  filter((condition_1 == 'Liver_Jansen2019_Hib' & condition_2 == 'Liver_Jansen2019_Act') | (condition_1=='Liver_PDExp_Fed' & condition_2 == 'Liver_PDExp_NotFed')) %>% 
  mutate(comparison = paste(condition_1,':',condition_2,sep = '')) %>% 
  group_by(gene_id) %>% 
  add_count() %>% 
  arrange(-n)

liver.overlap %>% filter(n==2)


## Muscle
extractConsequenceSummary(
  muscle.SwitchListAnalyzedSubset,
  asFractionTotal = FALSE,
  plotGenes=FALSE,
)
extractConsequenceEnrichment(muscle.SwitchListAnalyzedSubset) # Nothing significant 

extractSwitchOverlap(
  subsetSwitchAnalyzeRlist(muscle.SwitchListAnalyzedSubset,
                           str_detect(muscle.SwitchListAnalyzedSubset$isoformFeatures$condition_1,'Jansen2019_Act|PDExp_Fed') & 
                             str_detect(muscle.SwitchListAnalyzedSubset$isoformFeatures$condition_2,'Jansen2019_Hib|PDExp_NotFed')),
  filterForConsequences=TRUE,
  plotIsoforms = T
)

muscle.overlap <- muscle.PDExp.allSigSwitches %>% 
  filter(switchConsequencesGene==TRUE & gene_switch_q_value < 0.05) %>% 
  filter((condition_2 == 'Muscle_Jansen2019_Hib' & condition_1 == 'Muscle_Jansen2019_Act') | (condition_1=='Muscle_PDExp_Fed' & condition_2 == 'Muscle_PDExp_NotFed')) %>% 
  mutate(comparison = paste(condition_1,':',condition_2,sep = '')) %>% 
  group_by(gene_id) %>% 
  add_count() %>% 
  arrange(-n)

muscle.overlap %>% filter(n==2)

#ZC3H11A - different switch
switchPlot(muscle.SwitchListAnalyzedSubset, gene = 'gene-ZC3H11A',condition1 = 'Muscle_PDExp_Fed',condition2 = 'Muscle_PDExp_NotFed')
switchPlot(muscle.SwitchListAnalyzedSubset, gene = 'gene-ZC3H11A',condition1 = 'Muscle_Jansen2019_Hib',condition2 = 'Muscle_Jansen2019_Act')


#LDB3 - reversed switch!!
switchPlot(muscle.SwitchListAnalyzedSubset, gene = 'gene-LDB3',condition1 = 'Muscle_PDExp_Fed',condition2 = 'Muscle_PDExp_NotFed')
switchPlot(muscle.SwitchListAnalyzedSubset, gene = 'gene-LDB3',condition1 = 'Muscle_Jansen2019_Hib',condition2 = 'Muscle_Jansen2019_Act')

#PLS3 - Different switch
switchPlot(muscle.SwitchListAnalyzedSubset, gene = 'gene-PLS3',condition1 = 'Muscle_PDExp_Fed',condition2 = 'Muscle_PDExp_NotFed')
switchPlot(muscle.SwitchListAnalyzedSubset, gene = 'gene-PLS3',condition1 = 'Muscle_Jansen2019_Hib',condition2 = 'Muscle_Jansen2019_Act')



# Overlap with DE genes - any shared? -------------------------------------

adi.DEgenes <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/adipose_PWRes_Fed.vs.NotFed.tsv') %>% filter(IHW_pvalue < 0.05)
liv.DEgenes <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/liver_PWRes_Fed.vs.NotFed.tsv') %>% filter(IHW_pvalue < 0.05)
mus.DEgenes <- read_tsv('analysis/gene_level_rnaseq/pairwise_results/muscle_PWRes_Fed.vs.NotFed.tsv') %>% filter(IHW_pvalue < 0.05)


adipose.PDExp.allSigSwitches %>% filter(gene_id %in% adi.DEgenes$id) %>% filter(condition_1=='Fat_PDExp_Fed' & condition_2 == 'Fat_PDExp_NotFed')
liver.PDExp.allSigSwitches %>% filter(gene_id %in% liv.DEgenes$id) %>% filter(condition_1=='Liver_PDExp_Fed' & condition_2 == 'Liver_PDExp_NotFed')
muscle.PDExp.allSigSwitches %>% filter(gene_id %in% mus.DEgenes$id) %>% filter(condition_1=='Muscle_PDExp_Fed' & condition_2 == 'Muscle_PDExp_NotFed')

#



