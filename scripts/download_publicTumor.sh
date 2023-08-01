#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=2:00:00
#SBATCH --mem=16GB
#SBATCH -o /hpcfs/users/a1680844/20131906_HickeyT_JC_NormalBreast/%x_%j.out
#SBATCH -e /hpcfs/users/a1680844/20131906_HickeyT_JC_NormalBreast/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wenjun.liu@adelaide.edu.au


PROJ = /hpcfs/users/a1680844/20131906_HickeyT_JC_NormalBreast/data/raw/fastq

wget --no-passive-ftp ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR335/009/SRR3353249/SRR3353249.fastq.gz -O SRR3353249_GSM2112727_Tumor_P3_Vehicle_24_Hours_RNAseq_Homo_sapiens_RNA-Seq.fastq.gz
