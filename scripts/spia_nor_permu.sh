#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=12:00:00
#SBATCH --mem=16GB
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wenjun.liu@adelaide.edu.au

module load Anaconda3/5.1.0
source activate R_env

Rscript /hpcfs/users/a1680844/20131906_HickeyT_JC_NormalBreast/analysis/spia_nor_permu.R

conda deactivate

