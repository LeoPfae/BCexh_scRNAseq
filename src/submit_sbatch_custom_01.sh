#!/bin/bash
#
#SBATCH --job-name=BCexh_01
#SBATCH --partition=satpathy,owners,normal
#SBATCH --time=1-00:00:10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=21GB

ml R/4.2.0
ml hdf5

Rscript custom_01.R
