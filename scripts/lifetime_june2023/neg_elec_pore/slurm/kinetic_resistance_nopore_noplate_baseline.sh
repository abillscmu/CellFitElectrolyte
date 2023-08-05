#!/bin/bash
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -p cpu
#SBATCH -A venkvis
#SBATCH -J s-kin-b
#SBATCH -o neg_elec_logs/kinetic_resistance_nopore_noplate_baseline_%j.out
#SBATCH -e neg_elec_logs/kinetic_resistance_nopore_noplate_baseline_%j.err
#SBATCH -t 3-00:00
#SBATCH --mem 16G


julia --project=. scripts/lifetime_june2023/neg_elec_test_matrix/jl/kinetic_resistance_nopore_noplate.jl baseline

