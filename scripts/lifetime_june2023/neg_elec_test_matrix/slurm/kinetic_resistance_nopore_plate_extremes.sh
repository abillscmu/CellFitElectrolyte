#!/bin/bash
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -p cpu
#SBATCH -A venkvis
#SBATCH -J pl-kin-e
#SBATCH -o neg_elec_logs/kinetic_resistance_nopore_plate_extremes_%j.out
#SBATCH -e neg_elec_logs/kinetic_resistance_nopore_plate_extremes_%j.err
#SBATCH -t 3-00:00
#SBATCH --mem 16G


julia --project=. scripts/lifetime_june2023/neg_elec_test_matrix/jl/kinetic_resistance_nopore_plate.jl extremes

