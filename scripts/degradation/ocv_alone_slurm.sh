#!/bin/bash
#SBATCH -A venkvis
#SBATCH -p cpu
#SBATCH -n 1
#SBATCH -c 28
#SBATCH -t 1-00:00
#SBATCH --mem 32G
#SBATCH -J 🧢
#SBATCH --array=1,2,5,6,7,9,10,11,12,13,15,16,17,20,22,23,24,25,26,27,28,30
#SBATCH -o capacity%a.out
#SBATCH -e capacity%a.err

julia --project=. scripts/degradation/ocv_alone_slurm.jl