#!/bin/bash
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -p cpu
#SBATCH -A venkvis
#SBATCH -J pl-diff-a
#SBATCH -o neg_elec_logs/diffusive_resistance_nopore_plate_first_%j.out
#SBATCH -e neg_elec_logs/diffusive_resistance_nopore_plate_first_%j.err
#SBATCH -t 3-00:00
#SBATCH --mem 16G


julia --project=. scripts/lifetime_june2023/neg_elec_test_matrix/jl/diffusive_resistance_nopore_plate.jl first_cells

