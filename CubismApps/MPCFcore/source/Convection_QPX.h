/*
 *  *  Convection_QPX.h
 *   *  TryThis
 *    *
 *     *  Created by Costas Bekas on 1/24/13.
 *      *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *       *
 *        */

#pragma once

#ifndef _PRECDIV_
static const bool precdiv = false;
#else
static const bool precdiv = true;
#endif

/* #define _ACCURACY_DEBUG_ */ 

#ifdef _ACCURACY_DEBUG_
#include <cmath>
#include <limits>
#include <cstdio>
static const bool verbose = false;
static const double mytol = std::numeric_limits<Real>::epsilon() * 500;
#endif 

#include "check_errors.h"
#include "Convection_CPP.h"
#include "../../MPCFthread/source/Weno_CPP.h"
#include "../../MPCFthread/source/Weno_QPX.h"

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

class Convection_QPX : public Convection_CPP
{		
	__align(_ALIGNBYTES_) struct ScratchPad { Real tmp[TempSOA::NX][4];};

	void _sse_xweno_minus(float * const in, float * const out) const
	{
		enum { SX = InputSOA::PITCH };

		WenoQPX_MinusFunctor weno_ftor;
		
		for(int dy=0; dy<TempSOA::NY; ++dy)
		{
			float * const row = in + SX * dy;
			
			/* optionally we could skip the reading of W4 and C4 */
			for(int dx=0; dx<TempSOA::NX; dx+=4)
			{	
				const vector4double W4 = vec_ld(0 * sizeof(Real), row + dx);
				const vector4double C4 = vec_ld(4 * sizeof(Real), row + dx);
				const vector4double E4 = vec_ld(8 * sizeof(Real), row + dx);

				const vector4double a = vec_perm(W4, C4, vec_gpci(01234));
				const vector4double b = vec_perm(W4, C4, vec_gpci(02345));
				const vector4double c = vec_perm(W4, C4, vec_gpci(03456));
				const vector4double d = C4;
				const vector4double e = vec_perm(C4, E4, vec_gpci(01234));
				
				const vector4double reconstruction = weno_ftor._weno(a, b, c, d, e);
				
				vec_sta(reconstruction, 4 * (dx + TempSOA::PITCH * dy), out);
			}
		}
	}
	
	void _sse_xweno_plus(float * const in, float * const out) const
	{
		enum { SX = InputSOA::PITCH };
		
		WenoQPX_PlusFunctor weno_ftor;
		
		for(int dy=0; dy<TempSOA::NY; ++dy)
		{
			float * const row = in + SX * dy;
			
			/* optionally we could skip the reading of W4 and C4 */
			for(int dx=0; dx<TempSOA::NX; dx += 4)
			{	
				const vector4double W4 = vec_ld(0 * sizeof(Real), row + dx);
				const vector4double C4 = vec_ld(4 * sizeof(Real), row + dx);
				const vector4double E4 = vec_ld(8 * sizeof(Real), row + dx);
				
				const vector4double a = vec_perm(W4, C4, vec_gpci(02345));
				const vector4double b = vec_perm(W4, C4, vec_gpci(03456));
				const vector4double c = C4;
				const vector4double d = vec_perm(C4, E4, vec_gpci(01234));
				const vector4double e = vec_perm(C4, E4, vec_gpci(02345));
				
				const vector4double reconstruction = weno_ftor._weno(a, b, c, d, e);
				
				vec_sta(reconstruction, 4 * (dx + TempSOA::PITCH * dy), out);
			}
		}
	}
		
	template<typename Carrot, int start>
	void _sse_yweno(float * const in, float * const out) const
	{
		ScratchPad scratchpad;
		
		float * const tmp = (float *)&scratchpad.tmp[0][0];
		
		static const int SX = InputSOA::PITCH;
		
		Carrot weno_ftor;
		
		for(int dy=0; dy<TempSOA::NY; dy+=4)
		{
			float * const ptr = (start * SX) + &in[dy];
			
			for(int dx=0; dx<TempSOA::NX; ++dx)
			{			
				float * const entry = ptr + dx * SX;
				
				const vector4double a = vec_lda(0L, entry);
				const vector4double b = vec_lda(sizeof(Real) * SX, entry);
				const vector4double c = vec_lda(sizeof(Real) * 2 * SX, entry);
				const vector4double d = vec_lda(sizeof(Real) * 3 * SX, entry);
				const vector4double e = vec_lda(sizeof(Real) * 4 * SX, entry);
				
				const vector4double reconstruction = weno_ftor._weno(a, b, c, d, e);
				
				vec_sta(reconstruction, 0L, tmp + 4 * dx);
			}
			
			for(int dx=0; dx<TempSOA::NX-1; dx += 4)
			{
				float * const entry_in = tmp + 4 * dx;
				
				vector4double data0 = vec_lda(0L, entry_in);
				vector4double data1 = vec_lda(sizeof(Real) * 4, entry_in);
				vector4double data2 = vec_lda(sizeof(Real) * 4 * 2, entry_in);
				vector4double data3 = vec_lda(sizeof(Real) * 4 * 3, entry_in);
				
				_DIEGO_TRANSPOSE4(data0, data1, data2, data3);
				
				float * const entry_out = out + dx + TempSOA::PITCH * dy;
				
				vec_sta(data0, 0L, entry_out);
				vec_sta(data1, sizeof(Real) * TempSOA::PITCH, entry_out);
				vec_sta(data2, sizeof(Real) * TempSOA::PITCH * 2, entry_out);
				vec_sta(data3, sizeof(Real) * TempSOA::PITCH * 3, entry_out);
			}
			
			{
				enum { 
					A0 = TempSOA::NX-1,
					A1 = TempSOA::NX-1 + TempSOA::PITCH,
					A2 = TempSOA::NX-1 + TempSOA::PITCH * 2,
					A3 = TempSOA::NX-1 + TempSOA::PITCH * 3,
					B = TempSOA::PITCH
				};
				
				const int base = B * dy; 
				
				out[A0 + base] = scratchpad.tmp[TempSOA::NX-1][0];
				out[A1 + base] = scratchpad.tmp[TempSOA::NX-1][1];
				out[A2 + base] = scratchpad.tmp[TempSOA::NX-1][2];
				out[A3 + base] = scratchpad.tmp[TempSOA::NX-1][3];
			}
		}
	}
	
	template<typename Potato>
	void _sse_zweno(float * const a, float * const b, 
						  float * const c, float * const d, 
						  float * const e , float * const out) const
	{
		enum { 
			IN_STRIDE = InputSOA::PITCH,
			OUT_STRIDE = TempSOA::PITCH
		};
		
		WenoQPX<Potato> wenoqpx;
		
		for(int dy=0; dy<TempSOA::NY; dy++)
			wenoqpx.stride<TempSOA::NX - 1>(a + IN_STRIDE * dy, b + IN_STRIDE * dy, c + IN_STRIDE * dy, 
					d + IN_STRIDE * dy, e + IN_STRIDE * dy, out + OUT_STRIDE * dy);
	}
	
	void _xweno_minus(const InputSOA& in, TempSOA& out) 
	{ 
#ifdef _ACCURACY_DEBUG_
		{			
			InputSOA fake;			
			for(int iy=-3; iy < _BLOCKSIZE_ + 3; ++iy)
				for(int ix=-4; ix < _BLOCKSIZE_ + 4; ++ix)
					fake.ref(ix, iy) = drand48();
			
			TempSOA myrefoutput;
			
			Convection_CPP::_xweno_minus(fake, myrefoutput);
			_sse_xweno_minus(const_cast<float *>(fake.ptr(-4,0)), &out.ref(0,0)); 
			
			for(int iy=0; iy<TempSOA::NY; iy++)
				check_error(mytol, myrefoutput.ptr(0, iy), out.ptr(0, iy), TempSOA::NX);
			
			printf("checking wenox-minus fatto!\n"); 
		}
#endif
		_sse_xweno_minus(const_cast<float *>(in.ptr(-4,0)), &out.ref(0,0)); 
	}
	void _xweno_pluss(const InputSOA& in, TempSOA& out)
	{ 
#ifdef _ACCURACY_DEBUG_
		{
			InputSOA fake;
			
			for(int iy=-3; iy < _BLOCKSIZE_ + 3; ++iy)
				for(int ix=-4; ix < _BLOCKSIZE_ + 4; ++ix)
					fake.ref(ix, iy) = drand48();
			
			TempSOA myrefoutput;
			
			Convection_CPP::_xweno_pluss(fake, myrefoutput);
			_sse_xweno_plus(const_cast<float *>(fake.ptr(-4,0)), &out.ref(0,0)); 
			
			for(int iy=0; iy<TempSOA::NY; iy++)
				check_error(mytol, myrefoutput.ptr(0, iy), out.ptr(0, iy), TempSOA::NX);
			
			printf("checking wenox-plus fatto!\n"); 
		}
#endif
		_sse_xweno_plus(const_cast<float *>(in.ptr(-4,0)), &out.ref(0,0)); 
	}
	
	void _yweno_minus(const InputSOA& in, TempSOA& out)
	{ 
#ifdef _ACCURACY_DEBUG_
		{
			InputSOA fake;
			
			for(int iy=-3; iy < _BLOCKSIZE_ + 3; ++iy)
				for(int ix=-4; ix < _BLOCKSIZE_ + 3; ++ix)
					fake.ref(ix, iy) = drand48();
			
			TempSOA myrefoutput;
			Convection_CPP::_yweno_minus(fake, myrefoutput);
			_sse_yweno<WenoQPX_MinusFunctor, 0>(const_cast<float *>(fake.ptr(0,-3)), &out.ref(0,0)); 	
			
			for(int iy=0; iy<TempSOA::NY; iy++)
				check_error(mytol, myrefoutput.ptr(0, iy), out.ptr(0, iy), TempSOA::NX);
			
			printf("checking wenoy-minus fatto!\n"); 
		}
#endif
		_sse_yweno<WenoQPX_MinusFunctor, 0>(const_cast<float *>(in.ptr(0,-3)), &out.ref(0,0)); 
	}
	
	void _yweno_pluss(const InputSOA& in, TempSOA& out)
	{ 
#ifdef _ACCURACY_DEBUG_
		{
			InputSOA fake;
			
			for(int iy=-3; iy < _BLOCKSIZE_ + 3; ++iy)
				for(int ix=-4; ix < _BLOCKSIZE_ + 3; ++ix)
					fake.ref(ix, iy) = drand48();
			
			TempSOA myrefoutput;
			Convection_CPP::_yweno_pluss(fake, myrefoutput);
			_sse_yweno<WenoQPX_PlusFunctor, 1>(const_cast<float *>(fake.ptr(0,-3)), &out.ref(0,0)); 	
			
			for(int iy=0; iy<TempSOA::NY; iy++)
				check_error(mytol, myrefoutput.ptr(0, iy), out.ptr(0, iy), TempSOA::NX);
			
			printf("checking wenoy-plus fatto!\n"); 
		}
#endif
		_sse_yweno<WenoQPX_PlusFunctor, 1>(const_cast<float *>(in.ptr(0,-3)), &out.ref(0,0)); 
	} 
	
	void _zweno_minus(const int r, const RingInputSOA& in, TempSOA& out)
	{
#ifdef _ACCURACY_DEBUG_
		{
			RingInputSOA fake;
			
			for(int slice = 0 ; slice< 6 ; ++slice)
			for(int iy=-3; iy < _BLOCKSIZE_ + 3; ++iy)
				for(int ix=-4; ix < _BLOCKSIZE_ + 3; ++ix)
					fake.ref(slice).ref(ix, iy) = drand48();
			
			TempSOA myrefoutput;
			Convection_CPP::_zweno_minus(r, fake, myrefoutput);
			
			_sse_zweno<WenoQPX_MinusFunctor>(const_cast<float *>(fake(r-3).ptr(0,0)), 
											 const_cast<float *>(fake(r-2).ptr(0,0)), 
											 const_cast<float *>(fake(r-1).ptr(0,0)), 
											 const_cast<float *>(fake(r).ptr(0,0)), 
											 const_cast<float *>(fake(r+1).ptr(0,0)), 
											 &out.ref(0,0));
			
			for(int iy=0; iy<TempSOA::NY; iy++)
				check_error(mytol, myrefoutput.ptr(0, iy), out.ptr(0, iy), TempSOA::NX-1);
			
			printf("checking wenoz-minus fatto!\n");
		}
#endif
		_sse_zweno<WenoQPX_MinusFunctor>(const_cast<float *>(in(r-3).ptr(0,0)), 
						 const_cast<float *>(in(r-2).ptr(0,0)), 
						 const_cast<float *>(in(r-1).ptr(0,0)), 
						 const_cast<float *>(in(r).ptr(0,0)), 
						 const_cast<float *>(in(r+1).ptr(0,0)), 
						 &out.ref(0,0)); 
	}
	
	void _zweno_pluss(const int r, const RingInputSOA& in, TempSOA& out)
	 {
#ifdef _ACCURACY_DEBUG_
		 {
			 RingInputSOA fake;
			 
			 for(int slice = 0 ; slice< 6 ; ++slice)
				 for(int iy=-3; iy < _BLOCKSIZE_ + 3; ++iy)
					 for(int ix=-4; ix < _BLOCKSIZE_ + 3; ++ix)
						 fake.ref(slice).ref(ix, iy) = drand48();
			 
			 TempSOA myrefoutput;
			 Convection_CPP::_zweno_pluss(r, fake, myrefoutput);
			 
			 _sse_zweno<WenoQPX_PlusFunctor>(const_cast<float *>(fake(r-2).ptr(0,0)), 
											  const_cast<float *>(fake(r-1).ptr(0,0)), 
											  const_cast<float *>(fake(r).ptr(0,0)), 
											  const_cast<float *>(fake(r+1).ptr(0,0)), 
											  const_cast<float *>(fake(r+2).ptr(0,0)), 
											  &out.ref(0,0));
			 
			 for(int iy=0; iy<TempSOA::NY; iy++)
				 check_error(mytol, myrefoutput.ptr(0, iy), out.ptr(0, iy), TempSOA::NX-1);
			 
			 printf("checking wenoz-plus fatto!\n");
		 }
#endif
		 _sse_zweno<WenoQPX_PlusFunctor>(const_cast<float *>(in(r-2).ptr(0,0)), 
										 const_cast<float *>(in(r-1).ptr(0,0)), 
										 const_cast<float *>(in(r).ptr(0,0)), 
										 const_cast<float *>(in(r+1).ptr(0,0)), 
										 const_cast<float *>(in(r+2).ptr(0,0)), &out.ref(0,0)); 
	 } 
	
	void _sse_convert_aligned(float * const gptfirst, const int gptfloats, const int rowgpts,
							  float * const rho, float * const u, float * const v, float * const w, 
							  float * const p, float * const G, float * const P)
	{
		const vector4double F_1_2 = vec_splats(0.5f);
		const vector4double M_1_2 = vec_splats(-0.5f);
		const vector4double F_1 = vec_splats(1);
		
#define DESTID (dx + (InputSOA::PITCH)*dy)
		
		for(int dy=0; dy<_BLOCKSIZE_+6; dy++)
		{
			float * const in = gptfirst + dy*gptfloats*rowgpts -gptfloats;
			
			for(int dx=0; dx<_BLOCKSIZE_+8; dx+=4)
			{
				const int WID = (dx + (int)(dx==0))*gptfloats;
				const int CID = (dx+1)*gptfloats;
				const int EID = (dx+3 - (int)(dx==_BLOCKSIZE_+4))*gptfloats;
				
				vector4double dataA0 = vec_lda(0L, in + WID);
				vector4double dataA1 = vec_lda(0L, in + WID + 4);
				vector4double dataB0 = vec_lda(0L, in + CID );
				vector4double dataB1 = vec_lda(0L, in + CID + 4);
				vector4double dataC0 = vec_lda(0L, in + CID + gptfloats);
				vector4double dataC1 = vec_lda(0L, in + CID + gptfloats + 4);
				vector4double dataD0 = vec_lda(0L, in + EID );
				vector4double dataD1 = vec_lda(0L, in + EID + 4);
				
				_DIEGO_TRANSPOSE4(dataA0, dataB0, dataC0, dataD0);
				
				vec_sta(dataA0, 0L, rho + DESTID);
				
				const vector4double inv_rho = myreciprocal<precdiv>(dataA0);	

				vec_sta(vec_mul(dataB0, inv_rho), 0L, u + DESTID);
				vec_sta(vec_mul(dataC0, inv_rho), 0L, v + DESTID);
				vec_sta(vec_mul(dataD0, inv_rho), 0L, w + DESTID);
				
				_DIEGO_TRANSPOSE4(dataA1, dataB1, dataC1, dataD1);
				
				const vector4double alpha = vec_mul(M_1_2, inv_rho);
				const vector4double s_minus_P = vec_sub(dataA1, dataC1);
				const vector4double inv_G = myreciprocal<precdiv>(dataB1);	
				const vector4double speedsquared = vec_madd(dataB0, dataB0, vec_madd( dataC0, dataC0, vec_mul( dataD0, dataD0)));
				
				const vector4double myp = vec_mul(vec_madd(speedsquared, alpha, s_minus_P), inv_G);
				/*(dataA1 - (dataB0*dataB0 + dataC0*dataC0 + dataD0*dataD0)* (F_1_2*inv_rho) - dataC1)/dataB1);*/

				vec_sta(myp, 0L, p + DESTID);
				
				vec_sta(dataB1, 0L, G + DESTID);
				
				vec_sta(dataC1, 0L, P + DESTID);
			}
		}
		
#undef DESTID
		
	}
	
	/*void _convert(const Real * const gptfirst, const int gptfloats, const int rowgpts)
 * 	{
 * 			{
 * 						const size_t x = (size_t)gptfirst;
 * 									const int remainder = x & 0x1f;
 * 												
 * 															if (remainder)
 * 																		{
 * 																						printf("oooops! pointer is not aligned: 0x%x\n", (int)x);
 * 																										abort();
 * 																													}
 * 																															}
 * 																																	
 * 																																			_sse_convert_aligned(const_cast<Real*>(gptfirst), gptfloats, rowgpts, 
 * 																																										 & rho.ring.ref().ref(-4,-3),
 * 																																										 							 & u.ring.ref().ref(-4,-3),
 * 																																										 							 							 & v.ring.ref().ref(-4,-3),
 * 																																										 							 							 							 & w.ring.ref().ref(-4,-3),
 * 																																										 							 							 							 							 & p.ring.ref().ref(-4,-3),
 * 																																										 							 							 							 							 							 & G.ring.ref().ref(-4,-3), 
 * 																																										 							 							 							 							 							 							 & P.ring.ref().ref(-4,-3));
 * 																																										 							 							 							 							 							 							 				
 * 																																										 							 							 							 							 							 							 					}*/

public:
	Convection_QPX(const Real a, const Real dtinvh):
	Convection_CPP(a, dtinvh)
	{
	}
};
