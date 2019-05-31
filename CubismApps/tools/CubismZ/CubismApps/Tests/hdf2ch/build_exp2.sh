#!/bin/bash

set -x #echo on

if [ -z "$1" ]
then
	method=1
else
	method=$1
fi

bs=32
ds=512
nb=$(echo "$ds/$bs" | bc)

make clean

if   [ $method -eq 1 ]
then
	make all wavz=1 zlib=1
elif [ $method -eq 2 ]
then
	make all wavz=1 zlib=1 shuffle3=1
elif [ $method -eq 3 ]
then
	make all wavz=1 zlib=1 shuffle3=1 zerobits=4
elif [ $method -eq 4 ]
then
	make all wavz=1 zlib=1 shuffle3=1 zerobits=8
elif [ $method -eq 5 ]
then
	make all wavz=1 zlib=1 zshuffle=1
elif [ $method -eq 6 ]
then
	make all wavz=1 zlib=1 zshuffle=1 zerobits=4
elif [ $method -eq 7 ]
then
	make all wavz=1 zlib=1 zshuffle=1 zerobits=8
else
	echo "no valid option"
fi
