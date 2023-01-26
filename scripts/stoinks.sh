#!/bin/bash
#SBATCH -N 1
#SBATCH -J stoichiometry
#SBATCH -p RM-shared
#SBATCH -t 2-00:00
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem-per-cpu=2000M
#SBATCH --array=0-847

spack load julia@1.8
julia --project=. scripts/change_stoichiometry_ratio.jl
