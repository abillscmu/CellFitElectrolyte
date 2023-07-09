#!/bin/bash
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -p cpu
#SBATCH -A venkvis
#SBATCH -J 🎲🔋SEI 
#SBATCH -o thermal_sei.out
#SBATCH -e thermal_sei.err
#SBATCH -t 0-04:00
#SBATCH --mem 16G


julia --project=. scripts/lifetime_june2023/thermal_sei_buildup.jl

