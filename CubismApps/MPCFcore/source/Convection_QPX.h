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
				//(dataA1 - (dataB0*dataB0 + dataC0*dataC0 + dataD0*dataD0)* (F_1_2*inv_rho) - dataC1)/dataB1);

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

	inline void _qpx_xrhsadd(Real * const f, Real * const r)
	{
		const vector4double code = vec_gpci(01234);

		for(int iy=0; iy<OutputSOA::NY; iy++)
			for(int ix=0; ix<OutputSOA::NX; ix+=4)
			{
				Real * const entry_in = f + ix + TempSOA::PITCH*iy;

				const vector4double data0 = vec_lda(0L, entry_in);
				const vector4double data1 = vec_lda(sizeof(Real) * 4, entry_in);
				const vector4double datanext = vec_perm(data0, data1, code); 

				vec_sta(vec_sub(datanext, data0), 0L, r + ix + OutputSOA::PITCH*iy);
			}
	}

	void _xrhs()
	{	
		_qpx_xrhsadd(const_cast<Real*>(rho.flux().ptr(0,0)), &rho.rhs.ref(0,0));
		_qpx_xrhsadd(const_cast<Real*>(u.flux().ptr(0,0)), &u.rhs.ref(0,0));
		_qpx_xrhsadd(const_cast<Real*>(v.flux().ptr(0,0)), &v.rhs.ref(0,0));
		_qpx_xrhsadd(const_cast<Real*>(w.flux().ptr(0,0)), &w.rhs.ref(0,0));
		_qpx_xrhsadd(const_cast<Real*>(p.flux().ptr(0,0)), &p.rhs.ref(0,0));
		_qpx_xrhsadd(const_cast<Real*>(G.flux().ptr(0,0)), &G.rhs.ref(0,0));
		_qpx_xrhsadd(const_cast<Real*>(P.flux().ptr(0,0)), &P.rhs.ref(0,0));
	}

	inline void _qpx_yrhsadd(Real * const f, Real * const r)
	{
		enum { 
			SP = TempSOA::PITCH,
			DP = OutputSOA::PITCH
		} ;

		const vector4double code = vec_gpci(01234);

		for(int iy=0; iy<OutputSOA::NY; iy+=4)
			for(int ix=0; ix<OutputSOA::NX; ix+=4)
			{
				const int offset_in = iy + SP * ix;
				const int offset_out = ix + DP * iy;

				const vector4double dataA0 = vec_lda(0L, f + offset_in);
				const vector4double dataA1 = vec_lda(sizeof(Real) * 4, f + offset_in);
				const vector4double dataB0 = vec_lda(0L, f + offset_in + SP);
				const vector4double dataB1 = vec_lda(sizeof(Real) * 4, f + offset_in + SP);
				const vector4double dataC0 = vec_lda(0L, f + offset_in + 2 * SP);
				const vector4double dataC1 = vec_lda(sizeof(Real) * 4, f + offset_in + 2 * SP);
				const vector4double dataD0 = vec_lda(0L, f + offset_in + 3 * SP);
				const vector4double dataD1 = vec_lda(sizeof(Real) * 4, f + offset_in + 3 * SP);

				vector4double rhs0, rhs1, rhs2, rhs3;

				rhs0 = vec_sub(vec_perm(dataA0, dataA1, code), dataA0); 
				rhs1 = vec_sub(vec_perm(dataB0, dataB1, code), dataB0); 
				rhs2 = vec_sub(vec_perm(dataC0, dataC1, code), dataC0); 
				rhs3 = vec_sub(vec_perm(dataD0, dataD1, code), dataD0); 

				_DIEGO_TRANSPOSE4(rhs0, rhs1, rhs2, rhs3); 

				rhs0 = vec_add(rhs0, vec_lda(0L, r + offset_out));
				rhs1 = vec_add(rhs1, vec_lda(0L, r + offset_out + DP));
				rhs2 = vec_add(rhs2, vec_lda(0L, r + offset_out + 2 * DP));
				rhs3 = vec_add(rhs3, vec_lda(0L, r + offset_out + 3 * DP));

				vec_sta(rhs0, 0L, r + offset_out);
				vec_sta(rhs1, 0L, r + offset_out + DP);
				vec_sta(rhs2, 0L, r + offset_out + 2 * DP);
				vec_sta(rhs3, 0L, r + offset_out + 3 * DP);
			}
	}

	void _yrhs()
	{	
		_qpx_yrhsadd(const_cast<Real*>(rho.flux().ptr(0,0)), &rho.rhs.ref(0,0));
		_qpx_yrhsadd(const_cast<Real*>(u.flux().ptr(0,0)), &u.rhs.ref(0,0));
		_qpx_yrhsadd(const_cast<Real*>(v.flux().ptr(0,0)), &v.rhs.ref(0,0));
		_qpx_yrhsadd(const_cast<Real*>(w.flux().ptr(0,0)), &w.rhs.ref(0,0));
		_qpx_yrhsadd(const_cast<Real*>(p.flux().ptr(0,0)), &p.rhs.ref(0,0));
		_qpx_yrhsadd(const_cast<Real*>(G.flux().ptr(0,0)), &G.rhs.ref(0,0));
		_qpx_yrhsadd(const_cast<Real*>(P.flux().ptr(0,0)), &P.rhs.ref(0,0));
	}

	inline void _qpx_zrhsadd(Real * const fb, Real * const ff, Real * const r)
	{	
		for(int iy=0; iy<OutputSOA::NY; iy++)
			for(int ix=0; ix<OutputSOA::NX; ix+=4)
			{
				const vector4double result = vec_sub(vec_lda(0L, ff + ix + TempSOA::PITCH*iy), 
						vec_lda(0L, fb + ix + TempSOA::PITCH*iy));

				const vector4double oldvalue = vec_lda(0L, r + ix + OutputSOA::PITCH*iy);

				vec_sta(vec_add(result, oldvalue), 0L, r + ix + OutputSOA::PITCH*iy);
			}
	}

	void _zrhs()
	{	
		_qpx_zrhsadd(const_cast<Real*>(rho.flux(-1).ptr(0,0)), const_cast<Real*>(rho.flux(0).ptr(0,0)), &rho.rhs.ref(0,0));
		_qpx_zrhsadd(const_cast<Real*>(u.flux(-1).ptr(0,0)), const_cast<Real*>(u.flux(0).ptr(0,0)), &u.rhs.ref(0,0));
		_qpx_zrhsadd(const_cast<Real*>(v.flux(-1).ptr(0,0)), const_cast<Real*>(v.flux(0).ptr(0,0)), &v.rhs.ref(0,0));
		_qpx_zrhsadd(const_cast<Real*>(w.flux(-1).ptr(0,0)), const_cast<Real*>(w.flux(0).ptr(0,0)), &w.rhs.ref(0,0));
		_qpx_zrhsadd(const_cast<Real*>(p.flux(-1).ptr(0,0)), const_cast<Real*>(p.flux(0).ptr(0,0)), &p.rhs.ref(0,0));
		_qpx_zrhsadd(const_cast<Real*>(G.flux(-1).ptr(0,0)), const_cast<Real*>(G.flux(0).ptr(0,0)), &G.rhs.ref(0,0));
		_qpx_zrhsadd(const_cast<Real*>(P.flux(-1).ptr(0,0)), const_cast<Real*>(P.flux(0).ptr(0,0)), &P.rhs.ref(0,0));
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

	inline void _xsum(Real * const xm, Real * const xp, Real * const output)
	{
		for(int iy=0; iy<OutputSOA::NY; ++iy)
		{
			Real * const xmptr = xm + iy * TempSOA::PITCH;
			Real * const xpptr = xp + iy * TempSOA::PITCH;
			Real * const outputptr = output + iy * OutputSOA::PITCH;

			for(int ix=0; ix<OutputSOA::NX; ix += 4)
			{
				const vector4double datam = vec_perm(vec_lda(0L, xmptr + ix), 
						vec_lda(sizeof(Real) * 4, xmptr + ix), 
						vec_gpci(01234));				

				const vector4double datap = vec_lda(0L, xpptr + ix);

				vec_sta(vec_add(datam, datap), 0L, outputptr + ix);
			}
		}
	}

	inline vector4double _mixedterm(const vector4double am,  const vector4double ap,  
			const vector4double um,  const vector4double up)
	{
		const vector4double denom = myreciprocal<preclevel>(vec_sub(ap, am));
		const vector4double num = vec_msub(ap, um, vec_mul(am, up));

		return vec_mul(num, denom);
	}

	inline void _xsum(Real * const _am, Real * const _ap, Real * const _um, Real * const _up, Real * const output)
	{
		const vector4double code = vec_gpci(01234);

		for(int iy=0; iy<OutputSOA::NY; ++iy)
		{
			Real * const amptr = _am + iy * TempSOA::PITCH;
			Real * const apptr = _ap + iy * TempSOA::PITCH;
			Real * const umptr = _um + iy * TempSOA::PITCH;
			Real * const upptr = _up + iy * TempSOA::PITCH;

			Real * const outputptr = output + iy * OutputSOA::PITCH;

			for(int ix=0; ix<OutputSOA::NX; ix += 4)
			{
				const vector4double am0 = vec_lda(0L, amptr + ix);
				const vector4double am1 = vec_perm(am0, vec_lda(4 * sizeof(Real), amptr + ix), code);
				const vector4double ap0 = vec_lda(0L, apptr + ix);
				const vector4double ap1 = vec_perm(ap0, vec_lda(4 * sizeof(Real), apptr + ix), code);

				const vector4double um0 = vec_lda(0L, umptr + ix);
				const vector4double um1 = vec_perm(um0, vec_lda(4 * sizeof(Real), umptr + ix), code);
				const vector4double up0 = vec_lda(0L, upptr + ix);
				const vector4double up1 = vec_perm(up0, vec_lda(4 * sizeof(Real), upptr + ix), code);

				const vector4double result = vec_sub(_mixedterm(am1, ap1, um1, up1), _mixedterm(am0, ap0, um0, up0));

				vec_sta(result, 0L, outputptr + ix);
			}
		}
	}

	inline void _xextraterm(const TempSOA& um, const TempSOA& up, const TempSOA& Gm, const TempSOA& Gp
			, const TempSOA& Pm, const TempSOA& Pp
			, const TempSOA& am, const TempSOA& ap)
	{		
		_xsum(const_cast<Real *>(am.ptr(0,0)), const_cast<Real *>(ap.ptr(0,0)),
				const_cast<Real *>(um.ptr(0,0)), const_cast<Real *>(up.ptr(0,0)), &divu.ref(0,0));

		_xsum(const_cast<Real *>(Gm.ptr(0,0)), const_cast<Real *>(Gp.ptr(0,0)), &sumG.ref(0,0));
		_xsum(const_cast<Real *>(Pm.ptr(0,0)), const_cast<Real *>(Pp.ptr(0,0)), &sumP.ref(0,0));		
	}

	inline void _ysum(Real * const xm, Real * const xp, Real * const output)
	{
		enum { 
			SP = TempSOA::PITCH,
			DP = OutputSOA::PITCH,
			JUMP = sizeof(Real) * 4
		};

		const vector4double code = vec_gpci(01234);

		for(int iy=0; iy<OutputSOA::NY; iy+=4)
			for(int ix=0; ix<OutputSOA::NX; ix+=4)
			{
				const int offset_in0 = iy + SP * ix;
				const int offset_in1 = iy + SP * (ix + 1);
				const int offset_in2 = iy + SP * (ix + 2);
				const int offset_in3 = iy + SP * (ix + 3);

				const int offset_out = ix + DP * iy;

				vector4double rhs0 = vec_add(vec_ld(0L, xp + offset_in0), vec_perm(vec_lda(0L, xm + offset_in0), vec_lda(JUMP, xm + offset_in0), code)); 
				vector4double rhs1 = vec_add(vec_ld(0L, xp + offset_in1), vec_perm(vec_lda(0L, xm + offset_in1), vec_lda(JUMP, xm + offset_in1), code)); 
				vector4double rhs2 = vec_add(vec_ld(0L, xp + offset_in2), vec_perm(vec_lda(0L, xm + offset_in2), vec_lda(JUMP, xm + offset_in2), code)); 
				vector4double rhs3 = vec_add(vec_ld(0L, xp + offset_in3), vec_perm(vec_lda(0L, xm + offset_in3), vec_lda(JUMP, xm + offset_in3), code)); 

				_DIEGO_TRANSPOSE4(rhs0, rhs1, rhs2, rhs3); 

				rhs0 = vec_add(rhs0, vec_lda(0L, output + offset_out));
				rhs1 = vec_add(rhs1, vec_lda(0L, output + offset_out + DP));
				rhs2 = vec_add(rhs2, vec_lda(0L, output + offset_out + 2 * DP));
				rhs3 = vec_add(rhs3, vec_lda(0L, output + offset_out + 3 * DP));

				vec_sta(rhs0, 0L, output + offset_out);
				vec_sta(rhs1, 0L, output + offset_out + DP);
				vec_sta(rhs2, 0L, output + offset_out + 2 * DP);
				vec_sta(rhs3, 0L, output + offset_out + 3 * DP);
			}
	}

	inline void _ysum(Real * const am, Real * const ap,
			Real * const um, Real * const up, Real * const output)
	{
		enum { 
			SP = TempSOA::PITCH,
			DP = OutputSOA::PITCH,
			JUMP = sizeof(Real) * 4
		};

		TinyScratchPad scratchpad;

		Real * const tmp = (Real *)&scratchpad.tmp[0][0];

		const vector4double code = vec_gpci(01234);

		for(int iy=0; iy<OutputSOA::NY; iy+=4)
		{
			for(int ix=0; ix<OutputSOA::NX; ix+=4)
			{
				for(int j=0; j<4; ++j)
				{
					const int offset_in = iy + SP * (ix + j);

					const vector4double am0 = vec_ld(0L, am + offset_in);
					const vector4double am1 = vec_perm(am0, vec_ld(JUMP, am + offset_in), code);
					const vector4double ap0 = vec_ld(0L, ap + offset_in);
					const vector4double ap1 = vec_perm(ap0, vec_ld(JUMP, ap + offset_in), code);
					const vector4double um0 = vec_ld(0L, um + offset_in);
					const vector4double um1 = vec_perm(um0, vec_ld(JUMP, um + offset_in), code);
					const vector4double up0 = vec_ld(0L, up + offset_in);
					const vector4double up1 = vec_perm(up0, vec_ld(JUMP, up + offset_in), code);

					const vector4double result = vec_sub(_mixedterm(am1, ap1, um1, up1), _mixedterm(am0, ap0, um0, up0));

					vec_sta(result, 0L, tmp + 4 * j);
				}																					

				{
					vector4double sum0 = vec_lda(0L, tmp);
					vector4double sum1 = vec_lda(JUMP, tmp);
					vector4double sum2 = vec_lda(2 * JUMP, tmp);
					vector4double sum3 = vec_lda(3 * JUMP, tmp);

					_DIEGO_TRANSPOSE4(sum0, sum1, sum2, sum3); 

					const int offset_out = ix + DP * iy;

					sum0 = vec_add(sum0, vec_lda(0L, output + offset_out));
					sum1 = vec_add(sum1, vec_lda(0L, output + offset_out + DP));
					sum2 = vec_add(sum2, vec_lda(0L, output + offset_out + 2 * DP));
					sum3 = vec_add(sum3, vec_lda(0L, output + offset_out + 3 * DP));

					vec_sta(sum0, 0L, output + offset_out);
					vec_sta(sum1, 0L, output + offset_out + DP);
					vec_sta(sum2, 0L, output + offset_out + 2 * DP);
					vec_sta(sum3, 0L, output + offset_out + 3 * DP);
				}				
			}
		}
	}

	inline void _yextraterm(const TempSOA& um, const TempSOA& up, 
			const TempSOA& Gm, const TempSOA& Gp,
			const TempSOA& Pm, const TempSOA& Pp,
			const TempSOA& am, const TempSOA& ap)
	{
		_ysum(const_cast<Real *>(am.ptr(0,0)), const_cast<Real *>(ap.ptr(0,0)),
				const_cast<Real *>(um.ptr(0,0)), const_cast<Real *>(up.ptr(0,0)), &divu.ref(0,0));

		_ysum(const_cast<Real *>(Gm.ptr(0,0)), const_cast<Real *>(Gp.ptr(0,0)), &sumG.ref(0,0));
		_ysum(const_cast<Real *>(Pm.ptr(0,0)), const_cast<Real *>(Pp.ptr(0,0)), &sumP.ref(0,0));
	}

	inline void _zsum(Real * const xm, Real * const xp, Real * const output)
	{
		for(int iy=0; iy<OutputSOA::NY; ++iy)
		{
			Real * const xmptr = xm + iy * TempSOA::PITCH;
			Real * const xpptr = xp + iy * TempSOA::PITCH;
			Real * const outputptr = output + iy * OutputSOA::PITCH;

			for(int ix=0; ix<OutputSOA::NX; ix += 4)
			{
				const vector4double datam = vec_lda(0L, xmptr + ix);				
				const vector4double datap = vec_lda(0L, xpptr + ix);
				const vector4double oldvalue = vec_lda(0L, outputptr + ix);

				vec_sta(vec_add(oldvalue, vec_add(datam, datap)), 0L, outputptr + ix);
			}
		}
	}

	inline void _zsum(Real * const _am0, Real * const _ap0, Real * const _um0, Real * const _up0, 
			Real * const _am1, Real * const _ap1, Real * const _um1, Real * const _up1, 
			Real * const output)
	{
		for(int iy=0; iy<OutputSOA::NY; ++iy)
		{
			Real * const amptr0 = _am0 + iy * TempSOA::PITCH;
			Real * const apptr0 = _ap0 + iy * TempSOA::PITCH;
			Real * const umptr0 = _um0 + iy * TempSOA::PITCH;
			Real * const upptr0 = _up0 + iy * TempSOA::PITCH;

			Real * const amptr1 = _am1 + iy * TempSOA::PITCH;
			Real * const apptr1 = _ap1 + iy * TempSOA::PITCH;
			Real * const umptr1 = _um1 + iy * TempSOA::PITCH;
			Real * const upptr1 = _up1 + iy * TempSOA::PITCH;

			Real * const outputptr = output + iy * OutputSOA::PITCH;

			for(int ix=0; ix<OutputSOA::NX; ix += 4)
			{
				const vector4double am0 = vec_lda(0L, amptr0 + ix);
				const vector4double ap0 = vec_lda(0L, apptr0 + ix);
				const vector4double um0 = vec_lda(0L, umptr0 + ix);
				const vector4double up0 = vec_lda(0L, upptr0 + ix);

				const vector4double am1 = vec_lda(0L, amptr1 + ix);
				const vector4double ap1 = vec_lda(0L, apptr1 + ix);
				const vector4double um1 = vec_lda(0L, umptr1 + ix);
				const vector4double up1 = vec_lda(0L, upptr1 + ix);

				const vector4double result = vec_sub(_mixedterm(am1, ap1, um1, up1), _mixedterm(am0, ap0, um0, up0));
				const vector4double oldvalue = vec_lda(0L, outputptr + ix);

				vec_sta(vec_add(oldvalue, result), 0L, outputptr + ix);			
			}
		}
	}

	inline void _zextraterm(const TempSOA& um0, const TempSOA& up0, const TempSOA& um1, const TempSOA& up1, const TempSOA& Gm, const TempSOA& Gp
			, const TempSOA& Pm, const TempSOA& Pp
			, const TempSOA& am0, const TempSOA& ap0, const TempSOA& am1, const TempSOA& ap1)
	{
		_zsum(const_cast<Real *>(am0.ptr(0,0)), const_cast<Real *>(ap0.ptr(0,0)),
				const_cast<Real *>(um0.ptr(0,0)), const_cast<Real *>(up0.ptr(0,0)),
				const_cast<Real *>(am1.ptr(0,0)), const_cast<Real *>(ap1.ptr(0,0)),
				const_cast<Real *>(um1.ptr(0,0)), const_cast<Real *>(up1.ptr(0,0)),
				&divu.ref(0,0));

		_zsum(const_cast<Real *>(Gm.ptr(0,0)), const_cast<Real *>(Gp.ptr(0,0)), &sumG.ref(0,0));
		_zsum(const_cast<Real *>(Pm.ptr(0,0)), const_cast<Real *>(Pp.ptr(0,0)), &sumP.ref(0,0));
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
		{
			HLLESOA2D_QPX hllezator;

			hllezator.char_vel(rho.weno(0), rho.weno(1), u.weno(0), u.weno(1), p.weno(0), p.weno(1), G.weno(0), G.weno(1), P.weno(0), P.weno(1), charvel.ref(0), charvel.ref(1));
			hllezator.rho(rho.weno(0), rho.weno(1), u.weno(0), u.weno(1), charvel(0), charvel(1), rho.flux.ref());
			hllezator.pvel(rho.weno(0), rho.weno(1), u.weno(0), u.weno(1), p.weno(0), p.weno(1), charvel(0), charvel(1), u.flux.ref());
			hllezator.vel(rho.weno(0), rho.weno(1), v.weno(0), v.weno(1), u.weno(0), u.weno(1), charvel(0), charvel(1), v.flux.ref());
			hllezator.vel(rho.weno(0), rho.weno(1), w.weno(0), w.weno(1), u.weno(0), u.weno(1), charvel(0), charvel(1), w.flux.ref());     
			hllezator.e(rho.weno(0), rho.weno(1), u.weno(0), u.weno(1), v.weno(0), v.weno(1), w.weno(0), w.weno(1), p.weno(0), p.weno(1), G.weno(0), G.weno(1), P.weno(0), P.weno(1), charvel(0), charvel(1), p.flux.ref());       
			hllezator.rho(G.weno(0), G.weno(1), u.weno(0), u.weno(1), charvel(0), charvel(1), G.flux.ref());
			hllezator.rho(P.weno(0), P.weno(1), u.weno(0), u.weno(1), charvel(0), charvel(1), P.flux.ref());
		}

		_xextraterm(u.weno(0), u.weno(1), G.weno(0), G.weno(1), P.weno(0), P.weno(1), charvel(0), charvel(1));	
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

		{
			HLLESOA2D_QPX hllezator;

			hllezator.char_vel(rho.weno(0), rho.weno(1), v.weno(0), v.weno(1), p.weno(0), p.weno(1), G.weno(0), G.weno(1), P.weno(0), P.weno(1), charvel.ref(0), charvel.ref(1));
			hllezator.rho(rho.weno(0), rho.weno(1), v.weno(0), v.weno(1), charvel(0), charvel(1), rho.flux.ref());
			hllezator.vel(rho.weno(0), rho.weno(1), u.weno(0), u.weno(1), v.weno(0), v.weno(1), charvel(0), charvel(1), u.flux.ref());
			hllezator.pvel(rho.weno(0), rho.weno(1), v.weno(0), v.weno(1), p.weno(0), p.weno(1), charvel(0), charvel(1), v.flux.ref());
			hllezator.vel(rho.weno(0), rho.weno(1), w.weno(0), w.weno(1), v.weno(0), v.weno(1), charvel(0), charvel(1), w.flux.ref());     
			hllezator.e(rho.weno(0), rho.weno(1), v.weno(0), v.weno(1), u.weno(0), u.weno(1), w.weno(0), w.weno(1), p.weno(0), p.weno(1), G.weno(0), G.weno(1), P.weno(0), P.weno(1), charvel(0), charvel(1), p.flux.ref());
			hllezator.rho(G.weno(0), G.weno(1), v.weno(0), v.weno(1), charvel(0), charvel(1), G.flux.ref());       
			hllezator.rho(P.weno(0), P.weno(1), v.weno(0), v.weno(1), charvel(0), charvel(1), P.flux.ref());
		}

		_yextraterm(v.weno(0), v.weno(1), G.weno(0), G.weno(1), P.weno(0), P.weno(1), charvel(0), charvel(1));

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

		{
			HLLESOA2D_QPX hllezator;

			hllezator.char_vel(rho.weno(0), rho.weno(1), w.weno(0), w.weno(1), p.weno(0), p.weno(1), G.weno(0), G.weno(1), P.weno(0), P.weno(1), charvel.ref(0), charvel.ref(1));
			hllezator.rho(rho.weno(0), rho.weno(1), w.weno(0), w.weno(1), charvel(0), charvel(1), rho.flux.ref());
			hllezator.vel(rho.weno(0), rho.weno(1), u.weno(0), u.weno(1), w.weno(0), w.weno(1), charvel(0), charvel(1), u.flux.ref());
			hllezator.vel(rho.weno(0), rho.weno(1), v.weno(0), v.weno(1), w.weno(0), w.weno(1), charvel(0), charvel(1), v.flux.ref());
			hllezator.pvel(rho.weno(0), rho.weno(1), w.weno(0), w.weno(1), p.weno(0), p.weno(1), charvel(0), charvel(1), w.flux.ref());
			hllezator.e(rho.weno(0), rho.weno(1), w.weno(0), w.weno(1), u.weno(0), u.weno(1), v.weno(0), v.weno(1), p.weno(0), p.weno(1), G.weno(0), G.weno(1), P.weno(0), P.weno(1), charvel(0), charvel(1), p.flux.ref());
			hllezator.rho(G.weno(0), G.weno(1), w.weno(0), w.weno(1), charvel(0), charvel(1), G.flux.ref());    
			hllezator.rho(P.weno(0), P.weno(1), w.weno(0), w.weno(1), charvel(0), charvel(1), P.flux.ref());
		}

		_zextraterm(w.weno(-2), w.weno(-1),w.weno(0), w.weno(1), G.weno(0), G.weno(-1), P.weno(0), P.weno(-1), charvel(-2), charvel(-1), charvel(0), charvel(1));
	}

	public:
	Convection_QPX(const Real a, const Real dtinvh):
		Convection_CPP(a, dtinvh)
	{
	}
};
