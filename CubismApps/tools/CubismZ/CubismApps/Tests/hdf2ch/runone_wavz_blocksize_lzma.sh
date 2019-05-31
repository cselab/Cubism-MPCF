#!/bin/bash
set -x #echo on

##################
h5file=$1

#./genref_bs.sh $h5file 4
#make clean; make all wavz=1 lzma=1 shuffle3=1 bs=4
#./bench_wavz_blocksize.sh ./test_wavz_bs.sh $h5file 3 4

#./genref_bs.sh $h5file 8
#make clean; make all wavz=1 lzma=1 shuffle3=1 bs=8
#./bench_wavz_blocksize.sh ./test_wavz_bs.sh $h5file 3 8

./genref_bs.sh $h5file 16
make clean; make all wavz=1 lzma=1 shuffle3=1 bs=16
./bench_wavz_blocksize.sh ./test_wavz_bs.sh $h5file 3 16

./genref_bs.sh $h5file 32
make clean; make all wavz=1 lzma=1 shuffle3=1 bs=32
./bench_wavz_blocksize.sh ./test_wavz_bs.sh $h5file 3 32

./genref_bs.sh $h5file 64
make clean; make all wavz=1 lzma=1 shuffle3=1 bs=64
./bench_wavz_blocksize.sh ./test_wavz_bs.sh $h5file 3 64

exit 

./genref_bs.sh $h5file 128
make clean; make all wavz=1 lzma=1 shuffle3=1 bs=128
./bench_wavz_blocksize.sh ./test_wavz_bs.sh $h5file 3 128

./genref_bs.sh $h5file 256
make clean; make all wavz=1 lzma=1 shuffle3=1 bs=256
./bench_wavz_blocksize.sh ./test_wavz_bs.sh $h5file 3 256

./genref_bs.sh $h5file 512
make clean; make all wavz=1 lzma=1 shuffle3=1 bs=512
./bench_wavz_blocksize.sh ./test_wavz_bs.sh $h5file 3 512

