#!/bin/sh

#SBATCH --partition=gpu 
#SBATCH --gres=gpu:tesv100:1


ml conda
conda activate mNSF5

python fit_3samples_L10.py



