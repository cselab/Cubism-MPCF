/*
 * Author: Eric R. Schendel <erschend AT ncsu.edu>
 */

#include "isabela.h"

#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

static int isa_compress_buffer(char *in_buffer, int nx, int ny, int nz, double rate, int is_float, char *out_buffer, int *zbytes)
{
	enum ISABELA_status status;
	struct isabela_stream i_strm;
	struct isabela_config config;

	int element_bytes = sizeof(float);


	// Compress 1024 elements at a time
	config.window_size = 1024;

	// Approximate each window with 30 coefficients
	config.ncoefficients = 10;

	// Relative error between approximate and original values should be
	// no more than 5%
	config.error_rate = rate;

	// Size of each element
	config.element_byte_size = element_bytes;

	// Use either BSplines or Wavelets.
	 //config.transform = ISABELA_BSPLINES;
	config.transform = ISABELA_WAVELETS;

	int nbytes = nx*ny*nz*sizeof(float);

	// Input buffer

	// Setting output buffer with file_size
	//out_buffer = (char*) malloc (file_size);
	//assert (in_buffer != NULL && out_buffer != NULL);

	// Setup compression (deflate) with isabela_config
	status = isabelaDeflateInit (&i_strm, element_bytes, &config);
	assert (status == ISABELA_SUCCESS);

	i_strm.next_in = (void*) in_buffer;
	i_strm.avail_in = nbytes;
	i_strm.next_out = (void*) out_buffer;

	// Perform compression
	status = isabelaDeflate (&i_strm, ISABELA_FINISH);
	assert (status == ISABELA_SUCCESS);

//	printf ("\n");

//	printf ("%d bytes written out\n", i_strm.avail_out);
	static int ctr = 0;
//	printf ("%d bytes written out [%d]\n", i_strm.avail_out, ctr++);

	// Cleanup
	status = isabelaDeflateEnd (&i_strm);
	assert (status == ISABELA_SUCCESS);

	*zbytes = i_strm.avail_out;
	return 0;
}

static int isa_decompress_buffer(char *out_buffer, int nx, int ny, int nz, double rate, int is_float, char *in_buffer, int zbytes, int *uzbytes)
{
	enum ISABELA_status status;
	struct isabela_stream i_strm;
	struct isabela_config config;

	int element_bytes = sizeof(float);

	// Compress 1024 elements at a time
	config.window_size = 1024;

	// Approximate each window with 30 coefficients
	config.ncoefficients = 30;

	// Relative error between approximate and original values should be
	// no more than 5%
	config.error_rate = rate;

	// Size of each element
	config.element_byte_size = element_bytes;

	// Use either BSplines or Wavelets.
	config.transform = ISABELA_BSPLINES;
	// config.transform = ISABELA_WAVELETS;

	int nbytes = nx*ny*nz*sizeof(float);

	// Input buffer
	//in_buffer = (char*) malloc (file_size);

	// Being cautious and setting output buffer to 8 x compressed_file_size 
	//out_buffer = (char*) malloc (8 * file_size);
	//assert (in_buffer != NULL && out_buffer != NULL);

	// Read file
	//read_file (argv[1], in_buffer, file_size);

	// Setup deflation
	status = isabelaInflateInit (&i_strm, element_bytes, &config);
	assert (status == ISABELA_SUCCESS);

	// Setting up input / output buffers
	i_strm.next_in = (void*) in_buffer;
	i_strm.avail_in = zbytes;

	i_strm.next_out = (void*) out_buffer;

	status = isabelaInflate (&i_strm, ISABELA_FINISH);
	assert (status == ISABELA_SUCCESS);

//	printf ("\n");
//	printf ("%d bytes written out\n", i_strm.avail_out);
	*uzbytes = i_strm.avail_out;

	status = isabelaInflateEnd (&i_strm);
	assert (status == ISABELA_SUCCESS);

	static int ctr = 0;
	printf ("%d bytes written out [%d]\n", *uzbytes, ctr++);

	return 0;
}
