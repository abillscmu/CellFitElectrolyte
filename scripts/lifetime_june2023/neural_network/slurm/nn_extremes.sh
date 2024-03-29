#!/bin/bash
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -p cpu
#SBATCH -A venkvis
#SBATCH -J nn-e
#SBATCH -o neural_network_logs/extremes_%j.out
#SBATCH -e neural_network_logs/extremes_%j.err
#SBATCH -t 7-00:00
#SBATCH --mem 16G


julia --project=. scripts/lifetime_june2023/neural_network/jl/neural_network.jl extremes

