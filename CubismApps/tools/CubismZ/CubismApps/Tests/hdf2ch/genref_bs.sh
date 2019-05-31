#!/bin/bash
set -x #echo on

if [ -z "$1" ]
then
	echo "missing file"  
	exit
else
	h5file=$1
fi

if [ -z "$2" ]
then
	echo "missing bsize"  
	exit
else
	bsize=$2
fi

#bs=32
ds=512
nb=$(echo "$ds/$bsize" | bc)

rm -f ref.channel0

make clean
make all bs=$bsize

ulimit -s 25000000
export OMP_NUM_THREADS=1
./hdf2ch -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $h5file  -outdata c1 
mv c100000.StreamerGridPointIterative.channel0 ref.channel0




