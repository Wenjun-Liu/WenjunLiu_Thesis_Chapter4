#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=48:00:00
#SBATCH --mem=256GB
#SBATCH -o /hpcfs/users/a1680844/20131906_HickeyT_JC_NormalBreast/%x_%j.out
#SBATCH -e /hpcfs/users/a1680844/20131906_HickeyT_JC_NormalBreast/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wenjun.liu@adelaide.edu.au

## Cores
CORES=32

## Project Root
PROJ= /hpcfs/users/a1680844/20131906_HickeyT_JC_NormalBreast

## The environment containing snakemake
source activate Snakemake
cd /hpcfs/users/a1680844/20131906_HickeyT_JC_NormalBreast

## Run snakemake
snakemake \
  --cores ${CORES} \
  --use-conda \
  --notemp \
  --wrapper-prefix 'https://raw.githubusercontent.com/snakemake/snakemake-wrappers/'

## Add files to git
bash /hpcfs/users/a1680844/20131906_HickeyT_JC_NormalBreast/scripts/update_git.sh
