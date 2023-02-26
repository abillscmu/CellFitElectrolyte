#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=16000M
#SBATCH -J NUTS_3
#SBATCH -p RM-shared
#SBATCH -t 1-00:00
#SBATCH --array=0-847
#SBATCH -o logs/nuts%A_%a.out
#SBATCH -o logs/nuts%A_%a.err

/jet/home/abills/.julia/juliaup/bin/julia --project=. scripts/sandbox/newnuts_3.jl