bs ?= 16
vs ?= 128

MYFLAGS =  -D_BLOCKSIZE_=$(bs) -D_VOXELS_=$(vs) -D_FLOAT_PRECISION_ -D_SP_COMP_ -g -fopenmp -I../reader/ -I../../../Cubism/source/

vpcachier: main.o WaveletCompressor.o 
	$(CC) WaveletCompressor.o  main.o -o vpcachier  $(MYFLAGS) -lz -lstdc++

main.o: main.cpp
	$(CC) $(MYFLAGS) -c main.cpp -o main.o
	
WaveletCompressor.o: ../../MPCFnode/source/WaveletCompressor.cpp
	$(CC) $(MYFLAGS) -c ../../MPCFnode/source/WaveletCompressor.cpp -o WaveletCompressor.o
	
clean: *.o
	rm vpcachier *.o