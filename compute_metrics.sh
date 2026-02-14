#!/bin/bash
#SBATCH -c 1
#SBATCH --ntasks 3
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH -p cpu
#SBATCH -t 01:00:00
#SBATCH -o slurm-%j.out

make clean 
make 
./app
