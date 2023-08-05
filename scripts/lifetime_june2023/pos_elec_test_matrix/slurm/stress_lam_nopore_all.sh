#!/bin/bash
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -p cpu
#SBATCH -A venkvis
#SBATCH -J stress-all
#SBATCH -o pos_elec_logs/stress_all_files_%j.out
#SBATCH -e pos_elec_logs/stress_all_files_%j.err
#SBATCH -t 3-00:00
#SBATCH --mem 16G


julia --project=. scripts/lifetime_june2023/pos_elec_test_matrix/jl/stress_lam.jl all_files

