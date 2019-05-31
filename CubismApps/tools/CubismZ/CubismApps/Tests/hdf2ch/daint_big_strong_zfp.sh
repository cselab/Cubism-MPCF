#!/bin/bash -l
#

#SBATCH --job-name="cz_strong64"
#SBATCH --output=cubismz-%j.txt
#SBATCH --error=cubismz-%j.txt
#SBATCH --time=01:00:00
#SBATCH --nodes=64
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --constraint=gpu
#SBATCH --account=s659

#========================================
# load modules and run simulation
# export

export MPICH_GNI_DYNAMIC_CONN=disabled
export MPICH_GNI_USE_UNASSIGNED_CPUS=enabled
export MPICH_NEMESIS_ASYNC_PROGRESS=1
#export MPICH_MAX_THREAD_SAFETY=multiple

ulimit -s 1600000000
export OMP_NUM_THREADS=12
rm tmp00000.StreamerGridPointIterative.channel0; 

set -x #echo on

wt=3
inputfile=../../../../fabdata/data_010000-p-big.h5

thrval=0.005
rm -f cz_strong_zfp_5e3.txt
touch cz_strong_zfp_5e3.txt
srun --ntasks=4 --ntasks-per-node=1 ./hdf2ch -xpesize 1 -ypesize 2 -zpesize 2 -bpdx 64 -bpdy 32 -bpdz 32 -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt  >> cz_strong_zfp_5e3.txt
srun --ntasks=8 --ntasks-per-node=1 ./hdf2ch -xpesize 2 -ypesize 2 -zpesize 2 -bpdx 32 -bpdy 32 -bpdz 32 -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt  >> cz_strong_zfp_5e3.txt
srun --ntasks=16 --ntasks-per-node=1 ./hdf2ch -xpesize 2 -ypesize 2 -zpesize 4 -bpdx 32 -bpdy 32 -bpdz 16 -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> cz_strong_zfp_5e3.txt
srun --ntasks=32 --ntasks-per-node=1 ./hdf2ch -xpesize 2 -ypesize 4 -zpesize 4 -bpdx 32 -bpdy 16 -bpdz 16 -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> cz_strong_zfp_5e3.txt
srun --ntasks=64 --ntasks-per-node=1 ./hdf2ch -xpesize 4 -ypesize 4 -zpesize 4 -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> cz_strong_zfp_5e3.txt

thrval=0.1
rm -f cz_strong_zfp_1e1.txt
touch cz_strong_zfp_1e1.txt
srun --ntasks=4 --ntasks-per-node=1 ./hdf2ch -xpesize 1 -ypesize 2 -zpesize 2 -bpdx 64 -bpdy 32 -bpdz 32 -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt  >> cz_strong_zfp_1e1.txt
srun --ntasks=8 --ntasks-per-node=1 ./hdf2ch -xpesize 2 -ypesize 2 -zpesize 2 -bpdx 32 -bpdy 32 -bpdz 32 -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt  >> cz_strong_zfp_1e1.txt
srun --ntasks=16 --ntasks-per-node=1 ./hdf2ch -xpesize 2 -ypesize 2 -zpesize 4 -bpdx 32 -bpdy 32 -bpdz 16 -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> cz_strong_zfp_1e1.txt
srun --ntasks=32 --ntasks-per-node=1 ./hdf2ch -xpesize 2 -ypesize 4 -zpesize 4 -bpdx 32 -bpdy 16 -bpdz 16 -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> cz_strong_zfp_1e1.txt
srun --ntasks=64 --ntasks-per-node=1 ./hdf2ch -xpesize 4 -ypesize 4 -zpesize 4 -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> cz_strong_zfp_1e1.txt


thrval=0.6
rm -f cz_strong_zfp_6e1.txt
touch cz_strong_zfp_6e1.txt
srun --ntasks=4 --ntasks-per-node=1 ./hdf2ch -xpesize 1 -ypesize 2 -zpesize 2 -bpdx 64 -bpdy 32 -bpdz 32 -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt  >> cz_strong_zfp_6e1.txt
srun --ntasks=8 --ntasks-per-node=1 ./hdf2ch -xpesize 2 -ypesize 2 -zpesize 2 -bpdx 32 -bpdy 32 -bpdz 32 -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt  >> cz_strong_zfp_6e1.txt
srun --ntasks=16 --ntasks-per-node=1 ./hdf2ch -xpesize 2 -ypesize 2 -zpesize 4 -bpdx 32 -bpdy 32 -bpdz 16 -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> cz_strong_zfp_6e1.txt
srun --ntasks=32 --ntasks-per-node=1 ./hdf2ch -xpesize 2 -ypesize 4 -zpesize 4 -bpdx 32 -bpdy 16 -bpdz 16 -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> cz_strong_zfp_6e1.txt
srun --ntasks=64 --ntasks-per-node=1 ./hdf2ch -xpesize 4 -ypesize 4 -zpesize 4 -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt >> cz_strong_zfp_6e1.txt


