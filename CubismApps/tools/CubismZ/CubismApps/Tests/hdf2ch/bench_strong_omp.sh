#!/bin/bash -l
#

#SBATCH --job-name="cz_omp_strong"
#SBATCH --output=cz-%j.txt
#SBATCH --error=cz-%j.txt
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --constraint=gpu
#SBATCH --account=s659

#========================================
# load modules and run simulation
# export

#export MPICH_GNI_DYNAMIC_CONN=disabled
#export MPICH_GNI_USE_UNASSIGNED_CPUS=enabled
#export MPICH_NEMESIS_ASYNC_PROGRESS=1
#export MPICH_MAX_THREAD_SAFETY=multiple

rm -f cs_strong.txt
touch cs_strong.txt

set -x #echo on
#ulimit -s 1600000000
#export OMP_NUM_THREADS=1
rm tmp00000.StreamerGridPointIterative.channel0;

export CUBISMZ_NOIO=1
export OMP_PROC_BIND=TRUE
wt=3

inputfile=../../../../fabdata/data_010000-p.h5
nb=16

thrval=0.01
outfile=p10k_wt3_zlibB_t0.01.txt

export OMP_NUM_THREADS=1;  srun --ntasks=1 -c 1   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=2;  srun --ntasks=1 -c 2   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=4;  srun --ntasks=1 -c 4   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=6;  srun --ntasks=1 -c 6   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=8;  srun --ntasks=1 -c 8   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=10; srun --ntasks=1 -c 10  ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=12; srun --ntasks=1 -c 12  ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile

thrval=0.001
outfile=p10k_wt3_zlibB_t0.001.txt
export OMP_NUM_THREADS=1;  srun --ntasks=1 -c 1   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=2;  srun --ntasks=1 -c 2   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=4;  srun --ntasks=1 -c 4   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=6;  srun --ntasks=1 -c 6   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=8;  srun --ntasks=1 -c 8   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=10; srun --ntasks=1 -c 10  ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=12; srun --ntasks=1 -c 12  ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile

thrval=0.0001
outfile=p10k_wt3_zlibB_t0.0001.txt
export OMP_NUM_THREADS=1;  srun --ntasks=1 -c 1   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=2;  srun --ntasks=1 -c 2   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=4;  srun --ntasks=1 -c 4   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=6;  srun --ntasks=1 -c 6   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=8;  srun --ntasks=1 -c 8   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=10; srun --ntasks=1 -c 10  ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=12; srun --ntasks=1 -c 12  ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile


inputfile=../../../../fabdata/data_010000-p-med.h5
nb=32

thrval=0.01
outfile=p10k_wt3_zlibB_t0.01med.txt
export OMP_NUM_THREADS=1;  srun --ntasks=1 -c 1   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=2;  srun --ntasks=1 -c 2   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=4;  srun --ntasks=1 -c 4   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=6;  srun --ntasks=1 -c 6   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=8;  srun --ntasks=1 -c 8   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=10; srun --ntasks=1 -c 10  ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=12; srun --ntasks=1 -c 12  ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile

thrval=0.001
outfile=p10k_wt3_zlibB_t0.001med.txt
export OMP_NUM_THREADS=1;  srun --ntasks=1 -c 1   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=2;  srun --ntasks=1 -c 2   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=4;  srun --ntasks=1 -c 4   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=6;  srun --ntasks=1 -c 6   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=8;  srun --ntasks=1 -c 8   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=10; srun --ntasks=1 -c 10  ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=12; srun --ntasks=1 -c 12  ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile

thrval=0.0001
outfile=p10k_wt3_zlibB_t0.0001med.txt
export OMP_NUM_THREADS=1;  srun --ntasks=1 -c 1   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=2;  srun --ntasks=1 -c 2   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=4;  srun --ntasks=1 -c 4   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=6;  srun --ntasks=1 -c 6   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=8;  srun --ntasks=1 -c 8   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=10; srun --ntasks=1 -c 10  ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
export OMP_NUM_THREADS=12; srun --ntasks=1 -c 12  ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> $outfile
