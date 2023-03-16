#!/bin/bash

# Request 24 CPU cores
#SBATCH -N 1 -n 24
# Request maximum time
#SBATCH --time=7-00:00:00
#SBATCH -p proj_bhc
#SBATCH --mem 100G

echo $1
./kmonty 1e5 ../../../cks-run/dat-files/data0099.dat 1e19 $1
##4e18 for kappa runs
