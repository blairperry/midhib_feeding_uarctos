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
To add. 

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