!/bin/bash
NJOBS=$1

echo "#!/usr/bin/ksh
## submit with sbatch
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=06:30:00

#SBATCH --job-name="rosa-benchmarking-CVT"
#SBATCH --account=s70

export OMP_NUM_THREADS=32

aprun -n 1 -N 1 -d 32 python rosa-run-cvt.py gcc 390 /users/petrosk/projects/CubismApps/MPCFnode/benchmarks/cvt-ic/
" > submit-rosa-cvt

for I in $(seq 1 1 $NJOBS)
do 
sbatch submit-rosa-cvt
done 
