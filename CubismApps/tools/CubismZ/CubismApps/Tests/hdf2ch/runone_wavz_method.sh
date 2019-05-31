#!/bin/bash
set -x #echo on

##################
h5file=$1

./genref.sh $h5file

./build.sh 15   # make all wavz=1 zstd=1 shuffle3
./bench_wavz_method.sh ./test_wavz_type.sh $h5file 3  15

exit


./build.sh 0   # make all wavz=1 zlib=1 shuffle3
./bench_wavz_method.sh ./test_wavz_type.sh $h5file 3  0

./build.sh 14   # make all wavz=1 zstd=1 shuffle3
./bench_wavz_method.sh ./test_wavz_type.sh $h5file 3  14

exit


./build.sh 2   # make all wavz=1 lzma=1 shuffle3=1
./bench_wavz_method.sh ./test_wavz_type.sh $h5file 3 2

exit

./build.sh 0   # make all wavz=1 zlib=1 shuffle3
./bench_wavz_method.sh ./test_wavz_type.sh $h5file 3  0

./build.sh 1   # make all wavz=1 zlib=1
./bench_wavz_method.sh ./test_wavz_type.sh $h5file 3  1

./build.sh 2   # make all wavz=1 lzma=1 shuffle3=1
./bench_wavz_method.sh ./test_wavz_type.sh $h5file 3 2

./build.sh 3   # make all wavz=1 lzma=1 shuffle3=1 zerobits=4
./bench_wavz_method.sh ./test_wavz_type.sh $h5file 3 3

./build.sh 4   # make all wavz=1 lzma=1 shuffle3=1 zerobits=8
./bench_wavz_method.sh ./test_wavz_type.sh $h5file 3 4

./build.sh 5   # make all wavz=1 lzma=1 shuffle3=1 zerobits=12
./bench_wavz_method.sh ./test_wavz_type.sh $h5file 3 5

./build.sh 6   # make all wavz=1 lzma=1 shuffle3=1 zerobits=16
./bench_wavz_method.sh ./test_wavz_type.sh $h5file 3 6

./build.sh 13   # make all wavz=1 lzma=1
./bench_wavz_method.sh ./test_wavz_type.sh $h5file 3  13

