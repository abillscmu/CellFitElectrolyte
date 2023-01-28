#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -c 2
#SBATCH --mem=32000M
#SBATCH -J BNN-together
#SBATCH -p RM-shared
#SBATCH -t 2-00:00
#SBATCH -o logs/bnn_together%a.out
#SBATCH -e logs/bnn_together%a.err

stdbuf -i0 -o0 -e0 julia -p 16 --project=. scripts/degradation/degradation_new/BNN_thermal_simanneal_parallel_fullycoupled.jl