#!/bin/bash
set -x #echo on

if [ -z "$1" ]
then
    	echo "missing file"
        exit
else
    	h5file=$1
fi

if [ -z "$3" ]
then
	wt=3	
else
	wt=$3
fi

if [ -z "$4" ]
then
	bs=32	
else
	bs=$4
fi

#bs=32
ds=512
nb=$(echo "$ds/$bs" | bc)

rm tmp00000.StreamerGridPointIterative.channel0

ulimit -s 25000000
export OMP_NUM_THREADS=1
./hdf2ch -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $h5file -outdata tmp  -threshold $2 -wtype_write $wt

mpirun -n 8 ./ch2diff -simdata1 tmp00000.StreamerGridPointIterative.channel0  -simdata2 ref.channel0 -wtype $wt

