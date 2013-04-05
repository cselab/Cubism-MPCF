CC = mpic++
bs ?= 32

MYFLAGS = -O3 -mtune=native -march=native -ftree-vectorize -fstrict-aliasing -D_BLOCKSIZE_=$(bs) -D_FLOAT_PRECISION_ -D_SP_COMP_ -g -fopenmp -I../reader/ -I../../../Cubism/source

vp2hdf: main.o WaveletCompressor.o
	$(CC) WaveletCompressor.o  main.o  -o vp2hdf  $(MYFLAGS) -lhdf5 -lz -lstdc++

main.o: main.cpp
	$(CC) $(MYFLAGS) -c main.cpp -o main.o
	
WaveletCompressor.o: ../../MPCFnode/source/WaveletCompressor.cpp
	$(CC) $(MYFLAGS) -c ../../MPCFnode/source/WaveletCompressor.cpp -o WaveletCompressor.o
	
clean:
	rm vp2hdf *.o
