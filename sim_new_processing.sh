#!/bin/bash
#SBATCH --job-name=dbmt_process
#SBATCH --partition=fuchs
#SBATCH --nodes=3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=2000
#SBATCH --array=34-36
#SBATCH --mail-type=ALL
#SBATCH --time=00:05:00

srun Rscript sim_fit_process.R
exit 0
