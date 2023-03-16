#!/bin/bash
#SBATCH --job-name="bhac"
#SBATCH --output="bhac.%j.%N.out"
#SBATCH --partition=cca
##SBATCH --exclude=pcn-5-04,pcn-5-61,pcn-6-72
##SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --export=ALL
#SBATCH -t 2-00:00:00

#This job runs with 2 nodes, 16 cores per node for a total of 32 cores.
#ibrun in verbose mode will give binding detail

#module load intel/compiler/2017-4
#module load intel/mpi/2017-4
export I_MPI_EXTRA_FILESYSTEM=on
export I_MPI_EXTRA_FILESYSTEM_LIST=lustre
export FORT_BUFFERED=true

#####module load intel/license intel/compiler/2017-4 intel/mpi/2017-4 openmpi2/2.1.6-intel-hfi

for(( i=$1; i< $2; i=i+1))
do
	./kmonty 625000 ../raptor-cart/data0099.dat 1e19 $i
done

## -restart 170 -slice 171 -shell 171
