#!/bin/bash
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=16000M
#SBATCH -J stress-split
#SBATCH -p RM-shared
#SBATCH -t 2-00:00
#SBATCH --array=1,2,5,6,7,9,10,11,12,13,15,16,17,20,22,23,24,25,26,27,28,30
#SBATCH -o logs/stress_split%a.out
#SBATCH -e logs/stress_split%a.err

stdbuf -i0 -o0 -e0 julia --project=. scripts/degradation/degradation_new/Stress_alone_thermal_split.jl