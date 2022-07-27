#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -A venkvis
#SBATCH -p cpu
#SBATCH -t 2-00:00
#SBATCH -J OSGPP
#SBATCH --mem=12G
#SBATCH --array=1,2,5,6,7,9,10,11,12,13,15,16,17,20,22,23,24,25,26
#SBATCH --tmp=12G

/home/abills/julia-1.7.2/bin/julia reports/detailedreportbuilder.jl 
