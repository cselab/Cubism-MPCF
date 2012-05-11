#!/bin/bash
NJOBS=$1

echo "#!/bin/bash
## submit with sbatch
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=06:30:00

#SBATCH --job-name="todi-benchmarking-B2"

export OMP_NUM_THREADS=16

aprun -n 1 python todi-run.py gcc 390 128
" > submit-todi

for I in $(seq 1 1 $NJOBS)
do 
sbatch submit-todi
done
