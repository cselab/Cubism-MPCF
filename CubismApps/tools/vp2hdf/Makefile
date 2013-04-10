hdf = 1
include ../../Makefile.config 

CPPFLAGS += -I../reader 
CPPFLAGS += $(extra)

ifeq "$(zlib)" "1"
	MYFLAGS += -D_USE_ZLIB_
endif

ifeq "$(lz4)" "1"
	MYFLAGS += -D_USE_LZ4_  
endif    


vp2hdf: main.o WaveletCompressor.o
	$(CC) $(CPPFLAGS) WaveletCompressor.o  main.o  -o vp2hdf $(LIBS)

main.o: main.cpp
	$(CC) $(CPPFLAGS) -c main.cpp -o main.o
	
WaveletCompressor.o: ../../MPCFnode/source/WaveletCompressor.cpp
	$(CC) $(CPPFLAGS) -c ../../MPCFnode/source/WaveletCompressor.cpp -o WaveletCompressor.o
	
clean:
	rm -f vp2hdf *.o
