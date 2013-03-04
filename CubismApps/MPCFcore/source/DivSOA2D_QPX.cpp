#pragma once
#include "DivSOA2D_QPX.h"

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

inline void _qpx_xextraterm(Real * const xm, Real * const xp, Real * const output)
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

inline void _qpx_xextraterm(Real * const _am, Real * const _ap, Real * const _um, Real * const _up, Real * const output)
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


inline void _qpx_yextraterm(Real * const xm, Real * const xp, Real * const output)
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

__align(_ALIGNBYTES_) struct TinyScratchPadExtraTerm { Real tmp[4][4];};

inline void _qpx_yextraterm(Real * const am, Real * const ap,
		Real * const um, Real * const up, Real * const output)
{
	enum { 
		SP = TempSOA::PITCH,
		DP = OutputSOA::PITCH,
		JUMP = sizeof(Real) * 4
	};

	TinyScratchPadExtraTerm scratchpad;

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


inline void _qpx_zextraterm(Real * const xm, Real * const xp, Real * const output)
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

inline void _qpx_zextraterm(Real * const _am0, Real * const _ap0, Real * const _um0, Real * const _up0, 
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

void DivSOA2D_QPX::xrhs(const TempSOA& flux, OutputSOA& rhs) const
{
	_qpx_xrhsadd(const_cast<Real*>(flux.ptr(0,0)), &rhs.ref(0,0));
}

void DivSOA2D_QPX::yrhs(const TempSOA& flux, OutputSOA& rhs) const
{
	_qpx_yrhsadd(const_cast<Real*>(flux.ptr(0,0)), &rhs.ref(0,0));
}

void DivSOA2D_QPX::zrhs(const TempSOA& fback, const TempSOA& fforward, OutputSOA& rhs) const
{
	_qpx_zrhsadd(const_cast<Real*>(fback.ptr(0,0)), const_cast<Real*>(fforward.ptr(0,0)), &rhs.ref(0,0));
}

void DivSOA2D_QPX::xextraterm(const TempSOA& um, const TempSOA& up, const TempSOA& Gm, const TempSOA& Gp
		, const TempSOA& Pm, const TempSOA& Pp
		, const TempSOA& am, const TempSOA& ap,
		OutputSOA& divu, OutputSOA& sumG, OutputSOA& sumP)
{		
	_qpx_xextraterm(const_cast<Real *>(am.ptr(0,0)), const_cast<Real *>(ap.ptr(0,0)),
			const_cast<Real *>(um.ptr(0,0)), const_cast<Real *>(up.ptr(0,0)), &divu.ref(0,0));

	_qpx_xextraterm(const_cast<Real *>(Gm.ptr(0,0)), const_cast<Real *>(Gp.ptr(0,0)), &sumG.ref(0,0));
	_qpx_xextraterm(const_cast<Real *>(Pm.ptr(0,0)), const_cast<Real *>(Pp.ptr(0,0)), &sumP.ref(0,0));		
}

void DivSOA2D_QPX::yextraterm(const TempSOA& um, const TempSOA& up, 
		const TempSOA& Gm, const TempSOA& Gp,
		const TempSOA& Pm, const TempSOA& Pp,
		const TempSOA& am, const TempSOA& ap,
		OutputSOA& divu, OutputSOA& sumG, OutputSOA& sumP)
{
	_qpx_yextraterm(const_cast<Real *>(am.ptr(0,0)), const_cast<Real *>(ap.ptr(0,0)),
			const_cast<Real *>(um.ptr(0,0)), const_cast<Real *>(up.ptr(0,0)), &divu.ref(0,0));

	_qpx_yextraterm(const_cast<Real *>(Gm.ptr(0,0)), const_cast<Real *>(Gp.ptr(0,0)), &sumG.ref(0,0));
	_qpx_yextraterm(const_cast<Real *>(Pm.ptr(0,0)), const_cast<Real *>(Pp.ptr(0,0)), &sumP.ref(0,0));
}

void DivSOA2D_QPX::zextraterm(const TempSOA& um0, const TempSOA& up0, const TempSOA& um1, const TempSOA& up1, 
		const TempSOA& Gm, const TempSOA& Gp, const TempSOA& Pm, const TempSOA& Pp, 
		const TempSOA& am0, const TempSOA& ap0, const TempSOA& am1, const TempSOA& ap1,
		OutputSOA& divu, OutputSOA& sumG, OutputSOA& sumP)
{
	_qpx_zextraterm(const_cast<Real *>(am0.ptr(0,0)), const_cast<Real *>(ap0.ptr(0,0)),
			const_cast<Real *>(um0.ptr(0,0)), const_cast<Real *>(up0.ptr(0,0)),
			const_cast<Real *>(am1.ptr(0,0)), const_cast<Real *>(ap1.ptr(0,0)),
			const_cast<Real *>(um1.ptr(0,0)), const_cast<Real *>(up1.ptr(0,0)),
			&divu.ref(0,0));

	_qpx_zextraterm(const_cast<Real *>(Gm.ptr(0,0)), const_cast<Real *>(Gp.ptr(0,0)), &sumG.ref(0,0));
	_qpx_zextraterm(const_cast<Real *>(Pm.ptr(0,0)), const_cast<Real *>(Pp.ptr(0,0)), &sumP.ref(0,0));
}

