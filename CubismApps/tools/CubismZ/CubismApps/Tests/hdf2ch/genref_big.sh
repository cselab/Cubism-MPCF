#!/bin/bash -l
#

#SBATCH --job-name="genref64"
#SBATCH --output=cubismz-%j.txt
#SBATCH --error=cubismz-%j.txt
#SBATCH --time=00:10:00
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


#!/bin/bash
set -x #echo on

if [ -z "$1" ]
then
	echo "missing file"  
	exit
#h5file=../../../../fabdata/big/data_010000-p-big.h5
else
	h5file=$1
fi

#h5file=../../../../fabdata/big/data_010000-p-big.h5

bs=32
#ds=512
ds=2048
nb=$(echo "$ds/$bs" | bc)

#nb=2
rm -f ref.channel0

make clean
make all
export OMP_NUM_THREADS=1
#srun --ntasks=1 --ntasks-per-node=1 ./hdf2ch -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $h5file  -outdata c1 
srun --ntasks=64 --ntasks-per-node=1 ./hdf2ch -xpesize 4 -ypesize 4 -zpesize 4 -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata $h5file  -outdata c1 
mv c100000.StreamerGridPointIterative.channel0 ref.channel0




