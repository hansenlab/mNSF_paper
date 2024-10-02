#!/bin/sh

#SBATCH --partition=shared 

#SBATCH --cpus-per-task=1

ml conda
conda activate mNSF4
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/users/ywang/anaconda/envs/mNSF4/lib

python 02_model_fitting_disorted.py



