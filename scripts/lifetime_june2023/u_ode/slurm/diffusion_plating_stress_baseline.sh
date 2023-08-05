#!/bin/bash
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -p cpu
#SBATCH -A venkvis
#SBATCH -J u-base
#SBATCH -o uode_logs/baseline_%j.out
#SBATCH -e uode_logs/baseline_%j.err
#SBATCH -t 7-00:00
#SBATCH --mem 16G


julia --project=. scripts/lifetime_june2023/u_ode/jl/diffusion_plating_stress.jl baseline

