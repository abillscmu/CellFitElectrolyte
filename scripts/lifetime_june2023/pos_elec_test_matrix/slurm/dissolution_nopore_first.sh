#!/bin/bash
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -p cpu
#SBATCH -A venkvis
#SBATCH -J diss-first_cells
#SBATCH -o pos_elec_logs/diss_first_cells_%j.out
#SBATCH -e pos_elec_logs/diss_first_cells_%j.err
#SBATCH -t 0-24:00
#SBATCH --mem 16G


julia --project=. scripts/lifetime_june2023/pos_elec_test_matrix/jl/dissolution_nopore.jl first_cells

