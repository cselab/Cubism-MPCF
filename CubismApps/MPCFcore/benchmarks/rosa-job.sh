#!/bin/bash
NJOBS=$1

echo "#!/usr/bin/ksh
## submit with sbatch
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=00:390:00

#SBATCH --job-name="rosa-benchmarking-B1"
#SBATCH --reservation=main

export OMP_NUM_THREADS=32

aprun -n 1 -N 1 -d 32 python rosa-run.py gcc 390 128
" > submit-rosa

for I in $(seq 1 1 $NJOBS)
do 
sbatch submit-rosa
done 
