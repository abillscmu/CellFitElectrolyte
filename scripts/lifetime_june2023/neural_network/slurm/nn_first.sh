#!/bin/bash
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -p cpu
#SBATCH -A venkvis
#SBATCH -J nn-f
#SBATCH -o neural_network_logs/first_%j.out
#SBATCH -e neural_network_logs/first_%j.err
#SBATCH -t 7-00:00
#SBATCH --mem 16G


julia --project=. scripts/lifetime_june2023/neural_network/jl/neural_network.jl first_cells

