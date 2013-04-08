#pragma once

#ifdef _QPXEMU_
#include <xmmintrin.h>
#define vector4double __m128
#define _DIEGO_TRANSPOSE4(a,b,c,d) _MM_TRANSPOSE4_PS(a,b,c,d)
#define vec_lda(a, b) _mm_load_ps(b + (a)/4)
#define vec_sta(a, b, c) _mm_store_ps(c + (b)/4, a)
#define vec_mul(a, b) a * b
#define vec_sub(a, b) a - b
#define vec_madd(a, b, c) a * b + c
#define vec_nmsub(a, b, c) c - a * b
#define vec_splats(a) _mm_set1_ps(a)
#endif

#ifndef _DIEGO_TRANSPOSE4
#define _DIEGO_TRANSPOSE4(a, b, c, d)\
{\
const vector4double v01L = vec_perm(a, b, vec_gpci(00415));\
const vector4double v01H = vec_perm(a, b, vec_gpci(02637));\
const vector4double v23L = vec_perm(c, d, vec_gpci(00415));\
const vector4double v23H = vec_perm(c, d, vec_gpci(02637));\
\
a = vec_perm(v01L, v23L, vec_gpci(00145));\
b = vec_perm(v01L, v23L, vec_gpci(02367));\
c = vec_perm(v01H, v23H, vec_gpci(00145));\
d = vec_perm(v01H, v23H, vec_gpci(02367));\
}
#endif

namespace WaveletsOnInterval
{
struct WI4QPX
{
	template<int NH>
	struct __attribute__((__aligned__(_ALIGNBYTES_))) TinyScratchPad
	{
		FwtAp c[NH][4];
		FwtAp d[NH][4];
		
		void _cpbk(FwtAp *src, FwtAp * dst0, FwtAp * dst1, FwtAp * dst2, FwtAp * dst3)
		{
			for(int is = 0, id = 0; is < NH * 4; is += 16, id += 4)
			{
				vector4double d0 = vec_lda(0, src + is);
				vector4double d1 = vec_lda(sizeof(FwtAp) * 4, src + is);
				vector4double d2 = vec_lda(sizeof(FwtAp) * 4 * 2, src + is);
				vector4double d3 = vec_lda(sizeof(FwtAp) * 4 * 3, src + is);
				
				_DIEGO_TRANSPOSE4(d0, d1, d2, d3);
				
				vec_sta(d0, 0, dst0 + id);
				vec_sta(d1, 0, dst1 + id);
				vec_sta(d2, 0, dst2 + id);
				vec_sta(d3, 0, dst3 + id);
			}
		}
		
		void copyback(FwtAp * stream0, FwtAp * stream1, FwtAp * stream2, FwtAp * stream3)
		{
			assert(NH % 4 == 0);
			
			_cpbk(&c[0][0], stream0, stream1, stream2, stream3);
			_cpbk(&d[0][0], stream0 + NH, stream1 + NH, stream2 + NH, stream3 + NH);
		}
	};
	
	vector4double P_1_16;
	vector4double P_5_16;
	vector4double P_9_16;
	vector4double P_15_16;
	vector4double P_21_16;
	vector4double P_35_16;
	
	WI4QPX()
	{
		P_1_16 = vec_splats(1.f/16);
		P_5_16 = vec_splats(5.f/16);
		P_9_16 = vec_splats(9.f/16);
		P_15_16 = vec_splats(15.f/16);
		P_21_16 = vec_splats(21.f/16);
		P_35_16 = vec_splats(35.f/16);
	}
	inline vector4double interp_first(const vector4double f0, const vector4double f1, const vector4double f2, const vector4double f3)
	{
		return vec_madd(P_5_16, f0, vec_madd(P_15_16, f1, vec_nmsub(P_5_16, f2, vec_mul(P_1_16, f3))));
	};
	
	inline vector4double interp_middle(const vector4double f0, const vector4double f1, const vector4double f2, const vector4double f3)
	{
		return vec_nmsub(P_1_16, f0, vec_madd(P_9_16, f1, vec_nmsub(P_1_16, f3, vec_mul(P_9_16, f2))));
	};
	
	inline vector4double interp_onetolast(const vector4double f0, const vector4double f1, const vector4double f2, const vector4double f3)
	{
		return vec_madd(P_1_16, f0, vec_nmsub(P_5_16, f1, vec_madd(P_15_16, f2, vec_mul(P_5_16, f3))));
	};
	
	inline vector4double interp_last(const vector4double f0, const vector4double f1, const vector4double f2, const vector4double f3)
	{
		return vec_nmsub(P_5_16, f0, vec_madd(P_21_16, f1, vec_nmsub(P_35_16, f2, vec_mul(P_35_16, f3))));
	};
	
	template<int N>
	void forward(FwtAp * stream0, FwtAp * stream1, FwtAp * stream2, FwtAp * stream3)
	{		
		enum { NH = N/2 };

		assert(N >= 8);
		assert(N % 8 == 0);
		
		TinyScratchPad<NH> mycoeffs;
		
		vector4double d0 = vec_lda(0, stream0);
		vector4double d1 = vec_lda(0, stream1);
		vector4double d2 = vec_lda(0, stream2);
		vector4double d3 = vec_lda(0, stream3);
		
		_DIEGO_TRANSPOSE4(d0, d1, d2, d3);
		
		vector4double d4 = vec_lda(sizeof(FwtAp) * 4, stream0);
		vector4double d5 = vec_lda(sizeof(FwtAp) * 4, stream1);
		vector4double d6 = vec_lda(sizeof(FwtAp) * 4, stream2);
		vector4double d7 = vec_lda(sizeof(FwtAp) * 4, stream3);
		
		_DIEGO_TRANSPOSE4(d4, d5, d6, d7);
		
		vec_sta(d0, 0, mycoeffs.c[0]);
		vec_sta(vec_sub(d1, interp_first(d0, d2, d4, d6)), 0, mycoeffs.d[0]);
		vec_sta(d2, 0, mycoeffs.c[1]);
		vec_sta(vec_sub(d3, interp_middle(d0, d2, d4, d6)), 0, mycoeffs.d[1]);
		
		for(int src = 8, dst = 2; src < N; src += 4, dst += 2)
		{
			vector4double dm2 = d2;
			
			d0 = d4;
			d1 = d5;
			d2 = d6;
			d3 = d7;
			
			d4 = vec_lda(0, stream0 + src);
			d5 = vec_lda(0, stream1 + src);
			d6 = vec_lda(0, stream2 + src);
			d7 = vec_lda(0, stream3 + src);
			
			_DIEGO_TRANSPOSE4(d4, d5, d6, d7);
			
			vec_sta(d0, 0, mycoeffs.c[dst]);
			vec_sta(vec_sub(d1, interp_middle(dm2, d0, d2, d4)), 0, mycoeffs.d[dst]);
			vec_sta(d2, 0, mycoeffs.c[dst + 1]);
			vec_sta(vec_sub(d3, interp_middle(d0, d2, d4, d6)), 0, mycoeffs.d[dst + 1]);
		}
		
		vec_sta(d4, 0, mycoeffs.c[NH-2]);
		vec_sta(vec_sub(d5, interp_onetolast(d0, d2, d4, d6)), 0, mycoeffs.d[NH-2]);
		vec_sta(d6, 0, mycoeffs.c[NH-1]);
		vec_sta(vec_sub(d7, interp_last(d0, d2, d4, d6)), 0, mycoeffs.d[NH-1]);
		
		mycoeffs.copyback(stream0, stream1, stream2, stream3);
	}
};

template<int ROWSIZE, int COLSIZE>
struct WaveletSweepQPX
{
	WI4QPX qpxfwt;
	
	template<int BS, bool forward>
	inline void sweep1D(FwtAp data[BS][ROWSIZE])
	{
		if (forward)
			for(int iy = 0; iy < BS; iy += 4)
				qpxfwt.template forward<BS>(&data[iy][0], &data[iy + 1][0], &data[iy + 2][0], &data[iy + 3][0]);
		else 
			for(int iy = 0; iy < BS; ++iy)
				WI4<false>::template transform<BS, false>(&data[iy][0]);
	}
	
	template<int BS>
	inline void xy_transpose(FwtAp data[BS][ROWSIZE])
	{
		for(int iy = 0; iy < BS; ++iy)
			for(int ix = iy + 1; ix < BS; ++ix)
			{
				const FwtAp temp = data[iy][ix];
				data[iy][ix] = data[ix][iy];
				data[ix][iy] = temp;
			}
	}
	
	template<int BS>
	inline void xz_transpose(FwtAp data[BS][COLSIZE][ROWSIZE])
	{
		for(int iy = 0; iy < BS; ++iy)	
			for(int iz = 0; iz < BS; ++iz)
				for(int ix = iz + 1; ix < BS; ++ix)
				{
					const FwtAp temp = data[iz][iy][ix];
					data[iz][iy][ix] = data[ix][iy][iz];
					data[ix][iy][iz] = temp;
				}
	}
	
	template<int BS, bool forward>
	inline void sweep2D(FwtAp data[BS][ROWSIZE])
	{			
		sweep1D<BS, forward>(data);
		xy_transpose<BS>(data);
		sweep1D<BS, forward>(data);
	}
	
	template<int BS, bool bForward>
	inline void sweep3D(FwtAp data[BS][COLSIZE][ROWSIZE])
	{			
		if(bForward)
		{
			for(int iz = 0; iz < BS; ++iz)
				sweep2D<BS, true>(data[iz]);
			
			xz_transpose<BS>(data);
			
			for(int iz = 0; iz < BS; ++iz)
				for(int iy = 0; iy < BS; iy += 4)
					qpxfwt.template forward<BS>(&data[iz][iy][0], &data[iz][iy + 1][0], &data[iz][iy + 2][0], &data[iz][iy + 3][0]);
		}
		else
		{			
			for(int iz = 0; iz < BS; ++iz)
				for(int iy = 0; iy < BS; ++iy)
					WI4<false>::template transform<BS, false>(&data[iz][iy][0]);
			
			xz_transpose<BS>(data);
			
			for(int iz = 0; iz < BS; ++iz)
				sweep2D<BS, false>(data[iz]);
		}
	}
};
}

