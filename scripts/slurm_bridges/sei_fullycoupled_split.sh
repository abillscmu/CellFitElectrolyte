#!/bin/bash
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=16000M
#SBATCH -J SEI-split
#SBATCH -p RM-shared
#SBATCH -t 2-00:00
#SBATCH --array=1,2,5,6,7,9,10,11,12,13,15,16,17,20,22,23,24,25,26,27,28,30
#SBATCH -o logs/sei_split%a.out
#SBATCH -e logs/sei_split%a.err

julia --project=. scripts/degradation/degradation_new/SEI_alone_thermal_fullcouple_split.jl