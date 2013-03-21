/*
 *    Convection_QPX.h
 *     TryThis
 *    
 *       Created by Diego Rossinelli on 1/24/13.
 *        Copyright 2013 ETH Zurich. All rights reserved.
 *       
 */

#pragma once

#include "check_errors.h"
#include "Convection_CPP.h"
#include "WenoSOA2D_QPX.h"
#include "HLLESOA2D_QPX.h"
#include "DivSOA2D_QPX.h"

class Convection_QPX : public Convection_CPP
{		
	__align(_ALIGNBYTES_) struct TinyScratchPad { Real tmp[4][4];};
	
	void _qpx_convert_aligned(Real * const gptfirst, const int gptfloats, const int rowgpts,
							  Real * const rho, Real * const u, Real * const v, Real * const w, 
							  Real * const p, Real * const G, Real * const P)
	{
		const vector4double F_1_2 = vec_splats(0.5f);
		const vector4double M_1_2 = vec_splats(-0.5f);
		const vector4double F_1 = vec_splats(1);
		
#define DESTID (dx + (InputSOA::PITCH)*dy)
		
		for(int dy=0; dy<_BLOCKSIZE_+6; dy++)
		{
			Real * const in = gptfirst + dy*gptfloats*rowgpts -gptfloats;
			
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
				
				const vector4double inv_rho = myreciprocal<preclevel>(dataA0);	
				
				vec_sta(vec_mul(dataB0, inv_rho), 0L, u + DESTID);
				vec_sta(vec_mul(dataC0, inv_rho), 0L, v + DESTID);
				vec_sta(vec_mul(dataD0, inv_rho), 0L, w + DESTID);
				
				_DIEGO_TRANSPOSE4(dataA1, dataB1, dataC1, dataD1);
				
				const vector4double alpha = vec_mul(M_1_2, inv_rho);
				const vector4double s_minus_P = vec_sub(dataA1, dataC1);
				const vector4double inv_G = myreciprocal<preclevel>(dataB1);	
				const vector4double speedsquared = vec_madd(dataB0, dataB0, vec_madd( dataC0, dataC0, vec_mul( dataD0, dataD0)));
				
				const vector4double myp = vec_mul(vec_madd(speedsquared, alpha, s_minus_P), inv_G);
				
				vec_sta(myp, 0L, p + DESTID);
				vec_sta(dataB1, 0L, G + DESTID);
				vec_sta(dataC1, 0L, P + DESTID);
			}
		}
		
#undef DESTID
		
	}
	
	void _convert(const Real * const gptfirst, const int gptfloats, const int rowgpts)
	{
		{
			const size_t x = (size_t)gptfirst;
			const int remainder = x & 0x1f;
			assert(remainder == 0);
			assert(gptfloats == 16);
			if (remainder)
			{
				printf("oooops! pointer is not aligned: 0x%x\n", (int)x);
				abort();
			}
		}
		
		_qpx_convert_aligned(const_cast<Real*>(gptfirst), gptfloats, rowgpts, 
							 & rho.ring.ref().ref(-4,-3),
							 & u.ring.ref().ref(-4,-3),
							 & v.ring.ref().ref(-4,-3),
							 & w.ring.ref().ref(-4,-3),
							 & p.ring.ref().ref(-4,-3),
							 & G.ring.ref().ref(-4,-3), 
							 & P.ring.ref().ref(-4,-3));
		
	}
	
	void _xrhs()
	{
		DivSOA2D_QPX divtor;
		divtor.xrhs(rho.flux(), rho.rhs);
		divtor.xrhs(u.flux(), u.rhs);
		divtor.xrhs(v.flux(), v.rhs);
		divtor.xrhs(w.flux(), w.rhs);
		divtor.xrhs(p.flux(), p.rhs);
		divtor.xrhs(G.flux(), G.rhs);
		divtor.xrhs(P.flux(), P.rhs);	
	}
	
	void _yrhs()
	{
		DivSOA2D_QPX divtor;
		divtor.yrhs(rho.flux(), rho.rhs);
		divtor.yrhs(u.flux(), u.rhs);
		divtor.yrhs(v.flux(), v.rhs);
		divtor.yrhs(w.flux(), w.rhs);
		divtor.yrhs(p.flux(), p.rhs);
		divtor.yrhs(G.flux(), G.rhs);
		divtor.yrhs(P.flux(), P.rhs);	
	}
	
	void _zrhs()
	{
		DivSOA2D_QPX divtor;
		
		divtor.zrhs(rho.flux(-1), rho.flux(0), rho.rhs);
		divtor.zrhs(u.flux(-1), u.flux(0), u.rhs);
		divtor.zrhs(v.flux(-1), v.flux(0), v.rhs);
		divtor.zrhs(w.flux(-1), w.flux(0), w.rhs);
		divtor.zrhs(p.flux(-1), p.flux(0), p.rhs);
		divtor.zrhs(G.flux(-1), G.flux(0), G.rhs);
		divtor.zrhs(P.flux(-1), P.flux(0), P.rhs);	
	}
	
	void _copyback(Real * const gptfirst, const int gptfloats, const int rowgpts)
	{	
		const vector4double mya = vec_splats(a);
		const vector4double lambda = vec_splats(dtinvh);
		const vector4double M_1_6 = vec_splats(-1.f/6);
		
		const int offset1 = gptfloats;
		const int offset2 = 2 * gptfloats;
		const int offset3 = 3 * gptfloats;
		
		for(int iy=0; iy<OutputSOA::NY; iy++)
		{
			Real * const rhoptr = &rho.rhs.ref(0, iy);
			Real * const uptr = &u.rhs.ref(0, iy);
			Real * const vptr = &v.rhs.ref(0, iy);
			Real * const wptr = &w.rhs.ref(0, iy);
			Real * const pptr = &p.rhs.ref(0, iy);
			Real * const Gptr = &G.rhs.ref(0, iy);
			Real * const Pptr = &P.rhs.ref(0, iy);
			Real * const sumGptr = &sumG.ref(0, iy);
			Real * const sumPptr = &sumP.ref(0, iy);
			Real * const divuptr = &divu.ref(0, iy);
			
			for(int ix=0; ix<OutputSOA::NX; ix += 4)
			{
				const int entry_out = gptfloats*(ix + iy*rowgpts);
				
				vector4double d0 = vec_mul(lambda, vec_lda(0L, rhoptr + ix));
				vector4double d1 = vec_mul(lambda, vec_lda(0L, uptr + ix));
				vector4double d2 = vec_mul(lambda, vec_lda(0L, vptr + ix));
				vector4double d3 = vec_mul(lambda, vec_lda(0L, wptr + ix));
				vector4double d0B = vec_mul(lambda, vec_lda(0L, pptr + ix));
				vector4double d1B = vec_mul(M_1_6, vec_lda(0L, sumGptr + ix));
				vector4double d2B = vec_mul(M_1_6, vec_lda(0L, sumPptr + ix));
				vector4double d3B = vec_splats(0);
				
				vector4double mydivu = vec_lda(0L, divuptr + ix);
				d1B = vec_mul(lambda, vec_madd(d1, mydivu, vec_lda(0L, Gptr + ix)));
				d2B = vec_mul(lambda, vec_madd(d2, mydivu, vec_lda(0L, Pptr + ix)));
				
				_DIEGO_TRANSPOSE4(d0, d1, d2, d3);
				_DIEGO_TRANSPOSE4(d0B, d1B, d2B, d3B);
				
				d0	= vec_msub(mya, vec_lda(0L, gptfirst + entry_out), d0);
				d0B = vec_msub(mya, vec_lda(4 * sizeof(Real), gptfirst + entry_out), d0B);
				d1	= vec_msub(mya, vec_lda(0L, gptfirst + entry_out + offset1), d1);
				d1B = vec_msub(mya, vec_lda(4 * sizeof(Real), gptfirst + entry_out + offset1), d1B);
				d2	= vec_msub(mya, vec_lda(0L, gptfirst + entry_out + offset2), d2);
				d2B = vec_msub(mya, vec_lda(4 * sizeof(Real), gptfirst + entry_out + offset2), d2B);
				d3	= vec_msub(mya, vec_lda(0L, gptfirst + entry_out + offset3), d3);
				d3B	= vec_msub(mya, vec_lda(4 * sizeof(Real), gptfirst + entry_out + offset3), d3B);
				
				vec_sta(d0, 0L, gptfirst + entry_out);
				vec_sta(d0B, 4 * sizeof(Real), gptfirst + entry_out);
				vec_sta(d1, 0L, gptfirst + entry_out + offset1);
				vec_sta(d1B, 4 * sizeof(Real), gptfirst + entry_out + offset1);
				vec_sta(d2, 0L, gptfirst + entry_out + offset2);
				vec_sta(d2B, 4 * sizeof(Real), gptfirst + entry_out + offset2);
				vec_sta(d3, 0L, gptfirst + entry_out + offset3);					
				vec_sta(d3B, 4 * sizeof(Real), gptfirst + entry_out + offset3);
			}
		}
	}	
	
	void _xflux(const int relid)
	{
		{	
			WenoSOA2D_QPX wenoizer;
			
			wenoizer.xcompute(rho.ring(relid), rho.weno.ref(0), rho.weno.ref(1));
			wenoizer.xcompute(u.ring(relid), u.weno.ref(0), u.weno.ref(1));
			wenoizer.xcompute(v.ring(relid), v.weno.ref(0), v.weno.ref(1));
			wenoizer.xcompute(w.ring(relid), w.weno.ref(0), w.weno.ref(1));
			wenoizer.xcompute(p.ring(relid), p.weno.ref(0), p.weno.ref(1));
			wenoizer.xcompute(G.ring(relid), G.weno.ref(0), G.weno.ref(1));
			wenoizer.xcompute(P.ring(relid), P.weno.ref(0), P.weno.ref(1));
		}
		
		HLLESOA2D_QPX hllezator;
		hllezator.all(rho.weno(0), rho.weno(1), u.weno(0), u.weno(1), v.weno(0), v.weno(1), w.weno(0), w.weno(1), p.weno(0), p.weno(1), G.weno(0), G.weno(1), P.weno(0), P.weno(1), charvel.ref(0), charvel.ref(1), rho.flux.ref(), u.flux.ref(), v.flux.ref(), w.flux.ref(), p.flux.ref(), G.flux.ref(), P.flux.ref());
		
		DivSOA2D_QPX divtor;
		divtor.xextraterm(u.weno(0), u.weno(1), G.weno(0), G.weno(1), P.weno(0), P.weno(1), charvel(0), charvel(1), divu, sumG, sumP);	
	}
	
	void _yflux(const int relid)
	{	
		{
			WenoSOA2D_QPX wenoizer;
			
			wenoizer.ycompute(rho.ring(relid), rho.weno.ref(0), rho.weno.ref(1));
			wenoizer.ycompute(u.ring(relid), u.weno.ref(0), u.weno.ref(1));
			wenoizer.ycompute(v.ring(relid), v.weno.ref(0), v.weno.ref(1));
			wenoizer.ycompute(w.ring(relid), w.weno.ref(0), w.weno.ref(1));
			wenoizer.ycompute(p.ring(relid), p.weno.ref(0), p.weno.ref(1));
			wenoizer.ycompute(G.ring(relid), G.weno.ref(0), G.weno.ref(1));
			wenoizer.ycompute(P.ring(relid), P.weno.ref(0), P.weno.ref(1));
		}
		
		HLLESOA2D_QPX hllezator;
		hllezator.all(rho.weno(0), rho.weno(1), v.weno(0), v.weno(1), u.weno(0), u.weno(1), w.weno(0), w.weno(1), p.weno(0), p.weno(1), G.weno(0), G.weno(1), P.weno(0), P.weno(1), charvel.ref(0), charvel.ref(1), rho.flux.ref(), v.flux.ref(), u.flux.ref(), w.flux.ref(), p.flux.ref(), G.flux.ref(), P.flux.ref());
		
		DivSOA2D_QPX divtor;
		divtor.yextraterm(v.weno(0), v.weno(1), G.weno(0), G.weno(1), P.weno(0), P.weno(1), charvel(0), charvel(1), divu, sumG, sumP);
		
	}
	
	void _zflux(const int relid)
	{
		{
			WenoSOA2D_QPX wenoizer;
			
			wenoizer.zcompute(relid, rho.ring, rho.weno.ref(0), rho.weno.ref(1));
			wenoizer.zcompute(relid, u.ring, u.weno.ref(0), u.weno.ref(1));
			wenoizer.zcompute(relid, v.ring, v.weno.ref(0), v.weno.ref(1));
			wenoizer.zcompute(relid, w.ring, w.weno.ref(0), w.weno.ref(1));
			wenoizer.zcompute(relid, p.ring, p.weno.ref(0), p.weno.ref(1));
			wenoizer.zcompute(relid, G.ring, G.weno.ref(0), G.weno.ref(1));
			wenoizer.zcompute(relid, P.ring, P.weno.ref(0), P.weno.ref(1));	
		} 
		
		HLLESOA2D_QPX hllezator;
		hllezator.all(rho.weno(0), rho.weno(1), w.weno(0), w.weno(1), u.weno(0), u.weno(1), v.weno(0), v.weno(1), p.weno(0), p.weno(1), G.weno(0), G.weno(1), P.weno(0), P.weno(1), charvel.ref(0), charvel.ref(1), rho.flux.ref(), w.flux.ref(), u.flux.ref(), v.flux.ref(), p.flux.ref(), G.flux.ref(), P.flux.ref());		
		
		DivSOA2D_QPX divtor;
		divtor.zextraterm(w.weno(-2), w.weno(-1),w.weno(0), w.weno(1), G.weno(0), G.weno(-1), P.weno(0), P.weno(-1), charvel(-2), charvel(-1), charvel(0), charvel(1), divu, sumG, sumP);
	}
	
public:
	
	Convection_QPX(const Real a, const Real dtinvh): Convection_CPP(a, dtinvh) {}
};
