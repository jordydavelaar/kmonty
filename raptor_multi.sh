#!/bin/bash

# Request 24 CPU cores
#SBATCH -N 1 -n 24

# Request maximum time
#SBATCH --time=7-00:00:00

#SBATCH -p proj_bhc
#SBATCH --exclude=coma24
##SBATCH --mem 140G

./kmonty 51.2e8 ../../../library/sane_a+1o2/data1000.blk 4e19 0
