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


