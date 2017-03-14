#!/bin/bash

#SBATCH --job-name=neural_coding
#SBATCH --output=neural_coding.out
#SBATCH --error=neural_coding.err
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --account=mpcs51087

module load matlab/2016a
matlab neural_combos.m
