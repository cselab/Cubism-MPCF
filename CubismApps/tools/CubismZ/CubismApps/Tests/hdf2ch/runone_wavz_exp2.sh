#!/bin/bash
set -x #echo on

##################
h5file=$1

./genref.sh $h5file

#./build_exp2.sh 1   # make all wavz=1 zlib=1
#./bench_wavz_method.sh ./test_wavz_type.sh $h5file 3 1

#./build_exp2.sh 2   # make all wavz=1 zlib=1 shuffle3=1
#./bench_wavz_method.sh ./test_wavz_type.sh $h5file 3 2

#./build_exp2.sh 3   # make all wavz=1 zlib=1 shuffle3=1 zerobits=4
#./bench_wavz_method.sh ./test_wavz_type.sh $h5file 3 3

#./build_exp2.sh 4   # make all wavz=1 zlib=1 shuffle3=1 zerobits=8
#./bench_wavz_method.sh ./test_wavz_type.sh $h5file 3 4


./build_exp2.sh 5   # make all wavz=1 zlib=1 zshuffle=1
./bench_wavz_method.sh ./test_wavz_type.sh $h5file 3 2

./build_exp2.sh 6   # make all wavz=1 zlib=1 zshuffle=1 zerobits=4
./bench_wavz_method.sh ./test_wavz_type.sh $h5file 3 3

./build_exp2.sh 7   # make all wavz=1 zlib=1 zshuffle=1 zerobits=8
./bench_wavz_method.sh ./test_wavz_type.sh $h5file 3 4

