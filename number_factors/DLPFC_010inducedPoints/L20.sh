#!/bin/sh


#SBATCH --partition=gpu 
#SBATCH --gres=gpu:tesv100s:1

ml conda
conda activate mNSF5

python L20.py



