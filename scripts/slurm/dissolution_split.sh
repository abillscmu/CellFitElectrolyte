#!/bin/bash
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=16000M
#SBATCH -J diss-split
#SBATCH -p cpu
#SBATCH -t 2-00:00
#SBATCH --array=1,2,5,6,7,9,10,11,12,13,15,16,17,20,22,23,24,25,26,27,28,30
#SBATCH -o logs/diss_split%a.out
#SBATCH -e logs/diss_split%a.err

stdbuf -i0 -o0 -e0 julia --project=. scripts/degradation/degradation_new/dissolution_thermal_split.jl