#!/bin/bash
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -p cpu
#SBATCH -A venkvis
#SBATCH -J s-diff-s
#SBATCH -o neg_elec_logs/diffusive_resistance_nopore_noplate_first_cells_%j.out
#SBATCH -e neg_elec_logs/diffusive_resistance_nopore_noplate_first_cells_%j.err
#SBATCH -t 3-00:00
#SBATCH --mem 16G


julia --project=. scripts/lifetime_june2023/neg_elec_test_matrix/jl/diffusive_resistance_nopore_noplate.jl first_cells

