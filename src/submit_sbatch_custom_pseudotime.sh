#!/bin/bash
#
#SBATCH --job-name=BCexh_PT
#SBATCH --partition=satpathy
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=15GB

ml R/4.2.0
ml hdf5

Rscript pseudotime_reanalysis.R
