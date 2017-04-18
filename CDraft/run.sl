#!/bin/bash 
#SBATCH -n 1                      # Number of cores 
#SBATCH -N 1                      # Ensure that all cores are on one machine 
#SBATCH -t 4-0:00               # Runtime in D-HH:MM 
#SBATCH -p general                # Partition to submit to 
#SBATCH --mem=1000               # Memory pool for all cores (see also --mem-per-cpu) 
#SBATCH -o ./redfield.out
#SBATCH -e ./redfield.err
 
./redfield.x > result.out
