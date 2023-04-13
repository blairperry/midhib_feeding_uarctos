# Analysis of gene expression responses to mid-hibernation feeding in brown bears (*Ursus arctos*)

This repo contains code for all data processing and analysis used in our study of gene expression responses to mid-hibernation feeding in brown bears.

For questions, please contact blair.perry(at)wsu.edu.

## Data Availability
- Brown bear mRNA-seq data 
	- Feeding experiment data (generated in this study)
		- Raw data available at NCBI BioProject: [PRJNA640679](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA640679)
	- Hibernation versus active season data ([Jansen et al. 2019](https://www.nature.com/articles/s42003-019-0574-4))
		- Raw data available at NCBI BioProject: [PRJNA413091](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA413091)
	- Adipocyte culture experiment data ([Saxton et al. 2022](https://www.cell.com/iscience/fulltext/S2589-0042(22)01356-6))
		- Raw data are available at NCBI Bioproject [PRJNA578991](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA578991)

## Contents
[1. Quality trimming, mapping, and processing of RNA-seq data](#1-quality-trimming-mapping-and-processing-of-rna-seq-data)

IN PROGRESS

---

### 1. Quality trimming, mapping, and processing of RNA-seq data

Note: These analyses were run on WSU's HPC ([Kamiak](https://hpc.wsu.edu/)). For generalizability, simplified commands are presented here rather than the specific SLURM scripts used to run these commands on Kamiak. 

#### Trimming with TrimGalore
```bash
# Run TrimGalore
trim_galore --paired -q 20 --fastqc --fastqc_args "--noextract --nogroup --outdir 2_TrimGalore/fastqc/" --stringency 5 --illumina --length 50 -o trimmed_reads/ --clip_R1 12 --clip_R2 12 [path/to/read1] [path/to/read2]
```

#### Downsampling reads
- B4 sample has ~12M read pairs, while most others from that experiment have between 3-4M.
- Jansen 2019 samples have average of ~20-22M read pairs per sample.
- The feeding experiment samples (excluding B4) have on average 6.4M read pairs per sample. Going to downsample B4 and Jansen 2019 samples to that number. 

```bash
seqtk sample -s100 [path/to/trimmmed/read1] 6400000 > ./subsampled/[path/to/trimmed/read1/subset]
seqtk sample -s100 [path/to/trimmmed/read2] 6400000 > ./subsampled/[path/to/trimmed/read2/subset]
```

#### Mapping with STAR
Assembly: [GCF_023065955.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_023065955.1)

```bash
# Index genome for use with STAR
STAR --runMode genomeGenerate --runThreadN 8 --genomeDir ./star_reference --genomeFastaFiles GCF_023065955.1_UrsArc1.0_genomic.fna --sjdbGTFfile GCF_023065955.1_UrsArc1.0_genomic.gff

# Map Reads
STAR --genomeDir ./star_reference/ --runThreadN 8 --outFilterMultimapNmax 1 --twopassMode Basic --sjdbGTFfile GCF_023065955.1_UrsArc1.0_genomic.gff --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./star_mapped/[output_prefix] --readFilesIn [path/to/file1] [path/to/file2]
```

#### Quantifying gene-level read counts with featureCounts

```bash
# Convert gff to gtf w/ AGAT
agat_convert_sp_gff2gtf.pl --gff GCF_023065955.1_UrsArc1.0_genomic.gff -o GCF_023065955.1_UrsArc1.0_genomic.gtf

# Convert gtf to SAF using custom python script 
python2 ./utility_scripts/gtf_to_saf.allGeneInfo.py GCF_023065955.1_UrsArc1.0_genomic.gtf exon
# Output file is titled GCF_023065955.1_UrsArc1.0_genomic.allGeneInfo.saf

# Quantify gene-level counts using featureCounts
featureCounts -p -F 'gtf' -T 8 -t exon -g gene_id -a GCF_023065955.1_UrsArc1.0_genomic.allGeneInfo.saf -o ./postDexExperiment_08.10.22.txt [path/to/star_mapped]/*.sortedByCoord.out.bam
```


#### Transcipt-level quantification with Kallisto
Transcriptome: [GCF_023065955.1_UrsArc1.0_rna.fna](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/065/955/GCF_023065955.1_UrsArc1.0/GCF_023065955.1_UrsArc1.0_rna.fna.gz)

```bash
# Index transcriptome
kallisto index -i transcripts.idx GCF_023065955.1_UrsArc1.0_rna.fna

# Run Kallisto
kallisto quant -i transcripts.idx --rf-stranded -o kallisto_quants/$file_name [path/to/read1] [path/to/read2]
```


---

### 2. Gene-Level Expression Analyses

#### Gene-level Differential Expression Analyses (Tissue samples)

The following R script contains code used to:
-  Normalize gene expression counts for tissue-level samples of adipose, liver, and muscle (this study and Jansen et al. 2019)
-  Perform differential expression analyses between feeding experiment treatments and between active and hibernation season samples from Jansen et al. 2019

Link to Rscript: [_GeneLevel_DEseq2_08.10.22.R](https://github.com/blairperry/midhib_feeding_uarctos/blob/main/analysis/gene_level_rnaseq/_GeneLevel_DEseq2_08.10.22.R)

#### Gene Ontology (GO) and KEGG Pathway Overrepresentation Analysis of Differentially Expressed Genes

The following R script contains code used to:
- Characterization of enriched GO and KEGG terms for differentially expressed (DE) genes identified in above script

Link to Rscript: [DEGene_GOAnalysis_02.23.23.R](https://github.com/blairperry/midhib_feeding_uarctos/blob/main/analysis/gene_level_rnaseq/DEGene_GOAnalysis_02.23.23.R "DEGene_GOAnalysis_02.23.23.R")

#### Gene-level Differential Expression Analyses (Cell culture samples)

The following R script contains code used to:
-  Normalize gene expression counts for adipocyte culture samples (Saxton et al. 2022)
-  Perform differential expression analyses between HH (hibernation cells + hibernation serum) and HG (hibernation cells + post-feeding serum) treatments

Link to Rscript: [_CellCulture_GeneLevel_DEseq2_02.10.23.R](https://github.com/blairperry/midhib_feeding_uarctos/blob/main/analysis/gene_level_rnaseq/_CellCulture_GeneLevel_DEseq2_02.10.23.R "_CellCulture_GeneLevel_DEseq2_02.10.23.R")

#### Identification of genes with reversed expression after feeding

The following R script contains code used to:
-  Compare DE genes after mid-hibernation feeding with DE genes between active and hibernation seasons
-  Identify genes with "reversed" expression
	- Downregulated during hibernation -> upregulated after feeding
	- Upregulated during hibernation -> downregulated after feeding 
- Generate plot for Figure 1 of manuscript - venn diagrams showing overlap of DE genes from this study and Jansen et al 2019, dot-plots of reversed gene log2-fold changes
- GO and KEGG pathway enrichment analysis of reversed genes
	- Supplementary plotting of enrichment results

Link to Rscript: [2019vsPostDex_DEComparisons_11.01.22.R](https://github.com/blairperry/midhib_feeding_uarctos/blob/main/analysis/gene_level_rnaseq/2019vsPostDex_DEComparisons_11.01.22.R "2019vsPostDex_DEComparisons_11.01.22.R")

#### Comparison between tissue-level and cell culture gene expression 

The following R script contains code used to:
-  Compare differentially expressed (DE) genes after mid-hibernation feeding in adipose tissue with DE genes in hibernation adipocytes cultured with post-feeding serum (Saxton et al. 2022)
- Supplementary plotting of overlap and log2-fold changes of overlapping genes

Link to Rscript: [CellCultureVsPostDex_DEComparisons_02.13.23.R](https://github.com/blairperry/midhib_feeding_uarctos/blob/main/analysis/gene_level_rnaseq/CellCultureVsPostDex_DEComparisons_02.13.23.R "CellCultureVsPostDex_DEComparisons_02.13.23.R")

#### Prediction and Characterization of Upstream Regulatory Molecules with CHEA3

The following R script contains code used to:
-  Identify putative upstream reuglatory molecules for key subsets of DE genes and candidate serum proteins from Saxton et al 2022 using the CHEA3 API in R.

Link to Rscript: [chea3_API_Analysis_11.29.22.R](https://github.com/blairperry/midhib_feeding_uarctos/blob/main/analysis/gene_level_rnaseq/chea3_analyses/chea3_API_Analysis_11.29.22.R "chea3_API_Analysis_11.29.22.R")

The following R script contains code used to:
-  Identify the top predicted regulators for each gene set.
- Assess overlap of predicted regulators across tissues.
- Generate co-regulatory network plots for Figure 2 and supplementary figures
- Export network files for subsequent plotting in Cytoscape.

Link to Rscript: [chea3_resultWorkbook_11.29.22.R](https://github.com/blairperry/midhib_feeding_uarctos/blob/main/analysis/gene_level_rnaseq/chea3_analyses/chea3_resultWorkbook_11.29.22.R "chea3_resultWorkbook_11.29.22.R")