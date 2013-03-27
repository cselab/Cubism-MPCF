MYFLAGS =  -D_BLOCKSIZE_=32 -D_FLOAT_PRECISION_ -g

mytest: main.o WaveletCompressor.o
	$(CC) WaveletCompressor.o  main.o  -o mytest -lz -lstdc++

main.o: main.cpp
	$(CC) $(MYFLAGS) -c main.cpp -o main.o
	
WaveletCompressor.o: ../MPCFnode/source/WaveletCompressor.cpp
	$(CC) $(MYFLAGS) -c ../MPCFnode/source/WaveletCompressor.cpp -o WaveletCompressor.o