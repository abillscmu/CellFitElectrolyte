#!/bin/bash
#SBATCH -N 1
#SBATCH -n 15
#SBATCH --mem=0
#SBATCH -p idle
#SBATCH -A venkvis
#SBATCH -J 🧱🌊🔋
#SBATCH -t 7-00:00
#SBATCH -o logs/sei.out
#SBATCH -e logs/sei.err

julia --project=. -p 16 scripts/degradation/degradation_new/SEI_alone_thermal_fullcouple_parallel.jl