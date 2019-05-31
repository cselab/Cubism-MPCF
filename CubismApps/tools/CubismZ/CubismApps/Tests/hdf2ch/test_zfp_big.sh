#!/bin/bash
set -x #echo on

if [ -z "$1" ]
then
    	echo "missing file"
        exit
else
    	h5file=$1
fi

bs=32
ds=2048
nb=$(echo "$ds/$bs" | bc)
#wt=3

#make clean
#make all zfp=1

rm tmp00000.StreamerGridPointIterative.channel0

ulimit -s 160000000
export OMP_NUM_THREADS=1

srun --ntasks=64 --ntasks-per-node=1 ./hdf2ch -xpesize 4 -ypesize 4 -zpesize 4 -bpdx 16 -bpdy 16 -bpdz 16 -sim io -simdata $h5file -outdata tmp  -threshold $2
srun --ntasks=64 --ntasks-per-node=1 ./ch2diff -simdata1 tmp00000.StreamerGridPointIterative.channel0  -simdata2 ref.channel0


