#!/bin/sh

#SBATCH --partition=shared 

#SBATCH --cpus-per-task=5

ml conda_R/4.4

R CMD BATCH  step1_svg_eachSample.R

