#!/bin/bash 
#SBATCH --job-name=success
#SBATCH --time=00:30:00
#SBATCH --mail-type=end,fail
#SBATCH --mem=5g
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/lihanyu/slurm_output/slurm-%j.out
#SBATCH --account=wye2
#SBATCH --partition=standard
#SBATCH --array=1-500

module load Rtidyverse
R CMD BATCH --no-save --no-restore 6_600.R


