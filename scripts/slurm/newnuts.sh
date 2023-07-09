#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -p cpu
#SBATCH -A venkvis
#SBATCH -J 🥜
#SBATCH -t 1-00:00
#SBATCH --array=0-847
#SBATCH -o logs/mode%A_%a.out
#SBATCH -e logs/mode%A_%a.err

julia --project=. scripts/sandbox/newnuts.jl