#!/bin/bash

#SBATCH --job-name=neural_coding
#SBATCH --output=neural_coding.out
#SBATCH --error=neural_coding.err
#SBATCH --partition=broadwl
#SBATCH --ntasks=64
#SBATCH --account=mpcs51087

module load matlab
matlab par_neural_combos.m
