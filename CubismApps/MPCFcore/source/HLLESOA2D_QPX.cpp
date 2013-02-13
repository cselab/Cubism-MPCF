#include "HLLESOA2D_QPX.h"
#include "common.h"

template<int NX> inline void _qpx_hlle_rho(Real * const rm, Real * const rp,
		Real * const vm, Real * const vp,
		Real * const am, Real * const ap,
		Real * const out)
{
	enum { P = TempSOA::PITCH };

#define ID (ix + P*iy)	

	for(int iy=0; iy<TempSOA::NY; ++iy)
		for(int ix=0; ix<NX; ix+=4)
		{
			const vector4double aminus = vec_lda(0L, am + ID);
			const vector4double apluss = vec_lda(0L, ap + ID);

			const vector4double rminus = vec_lda(0L, rm + ID);
			const vector4double rpluss = vec_lda(0L, rp + ID);

			const vector4double fminus = vec_mul(vec_lda(0L, vm + ID), rminus);
			const vector4double fpluss = vec_mul(vec_lda(0L, vp + ID), rpluss);

			const vector4double tmp1 = vec_mul(vec_mul(aminus, apluss), vec_sub(rpluss, rminus));
			const vector4double tmp2 = vec_madd(apluss, fminus, vec_nmsub(aminus, fpluss, tmp1));
			const vector4double fother = vec_mul(tmp2, myreciprocal<preclevel>(vec_sub(apluss, aminus)));

			const vector4double result = vec_sel(fpluss, vec_sel(fother, fminus, aminus), apluss);

			vec_sta(result, 0L, out + ID);
		}
#undef ID
}

template<int NX> inline void _qpx_hlle_vel(Real * const rm, Real * const rp,
		Real * const vm, Real * const vp,
		Real * const vdm, Real * const vdp,
		Real * const am, Real * const ap,
		Real * const out)
{
	enum { P = TempSOA::PITCH };

#define ID (ix + P*iy)	

	for(int iy=0; iy<TempSOA::NY; ++iy)
		for(int ix=0; ix < NX; ix+=4)
		{
			const vector4double aminus = vec_lda(0L, am + ID);
			const vector4double apluss = vec_lda(0L, ap + ID);

			const vector4double uminus = vec_mul(vec_lda(0L, vm + ID), vec_lda(0L, rm + ID));
			const vector4double upluss = vec_mul(vec_lda(0L, vp + ID), vec_lda(0L, rp + ID));

			const vector4double fminus = vec_mul(vec_lda(0L, vdm + ID), uminus);
			const vector4double fpluss = vec_mul(vec_lda(0L, vdp + ID), upluss);

			const vector4double tmp1 = vec_mul(vec_mul(aminus, apluss), vec_sub(upluss, uminus));
			const vector4double tmp2 = vec_madd(apluss, fminus, vec_nmsub(aminus, fpluss, tmp1));
			const vector4double fother = vec_mul(tmp2, myreciprocal<preclevel>(vec_sub(apluss, aminus)));

			const vector4double result = vec_sel(fpluss, vec_sel(fother, fminus, aminus), apluss);

			vec_sta(result, 0L, out + ID);
		}
#undef ID
}


template<int NX> inline void _qpx_hlle_pvel(Real * const rm, Real * const rp,
		Real * const vm, Real * const vp,
		Real * const pm, Real * const pp,
		Real * const am, Real * const ap,
		Real * const out)
{
	enum { P = TempSOA::PITCH };

#define ID (ix + P*iy)
	for(int iy=0; iy<TempSOA::NY; ++iy)
		for(int ix=0; ix<NX; ix+=4)
		{			
			const vector4double myvminus = vec_lda(0L, vm + ID);
			const vector4double myvpluss = vec_lda(0L, vp + ID);

			const vector4double uminus = vec_mul(myvminus, vec_lda(0L, rm + ID));
			const vector4double upluss = vec_mul(myvpluss, vec_lda(0L, rp + ID));

			const vector4double fminus = vec_madd(myvminus, uminus, vec_lda(0L, pm + ID));
			const vector4double fpluss = vec_madd(myvpluss, upluss, vec_lda(0L, pp + ID));

			const vector4double aminus = vec_lda(0L, am + ID);
			const vector4double apluss = vec_lda(0L, ap + ID);

			const vector4double tmp1 = vec_mul(vec_mul(aminus, apluss), vec_sub(upluss, uminus));
			const vector4double tmp2 = vec_madd(apluss, fminus, vec_nmsub(aminus, fpluss, tmp1));
			const vector4double fother = vec_mul(tmp2, myreciprocal<preclevel>(vec_sub(apluss, aminus)));

			const vector4double result = vec_sel(fpluss, vec_sel(fother, fminus, aminus), apluss);

			vec_sta(result, 0L, out + ID);
		}
#undef ID
}

template<int NX> inline void _qpx_hlle_e(Real * const rm, Real * const rp,
		Real * const vdm, Real * const vdp,
		Real * const v1m, Real * const v1p,
		Real * const v2m, Real * const v2p,
		Real * const pm, Real * const pp,
		Real * const Gm, Real * const Gp, 
		Real * const PIm, Real * const PIp,
		Real * const am, Real * const ap,
		Real * const out)
{
	enum { P = TempSOA::PITCH };
#define ID (ix + P*iy)

	const vector4double F_1_2 = vec_splats(0.5);

	for(int iy=0; iy<TempSOA::NY; ++iy)
		for(int ix=0; ix<NX; ix+=4)
		{			
			const vector4double rminus	= vec_lda(0L, rm + ID);
			const vector4double vdminus	= vec_lda(0L, vdm + ID);
			const vector4double v1minus	= vec_lda(0L, v1m + ID);
			const vector4double v2minus	= vec_lda(0L, v2m + ID);
			const vector4double pminus	= vec_lda(0L, pm + ID);
			const vector4double Gminus	= vec_lda(0L, Gm + ID);
			const vector4double PIminus	= vec_lda(0L, PIm + ID);

			const vector4double speedminus = vec_madd(vdminus, vdminus, vec_madd(v1minus, v1minus, vec_mul(v2minus, v2minus)));
			const vector4double eminus = vec_madd(pminus, Gminus, vec_madd(vec_mul(F_1_2, rminus), speedminus, PIminus));

			const vector4double rplus	= vec_lda(0L, rp + ID);
			const vector4double vdplus	= vec_lda(0L, vdp + ID);
			const vector4double v1plus	= vec_lda(0L, v1p + ID);
			const vector4double v2plus	= vec_lda(0L, v2p + ID);
			const vector4double pplus	= vec_lda(0L, pp + ID);
			const vector4double Gplus	= vec_lda(0L, Gp + ID);
			const vector4double PIplus	= vec_lda(0L, PIp + ID);

			const vector4double speedplus = vec_madd(vdplus, vdplus, vec_madd(v1plus, v1plus, vec_mul(v2plus, v2plus)));
			const vector4double eplus = vec_madd(pplus, Gplus, vec_madd(vec_mul(F_1_2, rplus), speedplus, PIplus));

			const vector4double fminus = vec_mul(vdminus, vec_add(pminus, eminus));
			const vector4double fpluss = vec_mul(vdplus, vec_add(pplus, eplus));

			const vector4double aminus	= vec_lda(0L, am + ID);
			const vector4double apluss	= vec_lda(0L, ap + ID);

			const vector4double tmp1 = vec_mul(vec_mul(aminus, apluss), vec_sub(eplus, eminus));
			const vector4double tmp2 = vec_madd(apluss, fminus, vec_nmsub(aminus, fpluss, tmp1));
			const vector4double fother = vec_mul(tmp2, myreciprocal<preclevel>(vec_sub(apluss, aminus)));

			const vector4double result = vec_sel(fpluss, vec_sel(fother, fminus, aminus), apluss);

			vec_sta(result, 0L, out + ID);
		}
#undef ID
}

template<int NX> inline void _qpx_char_vel(Real * const rm, Real * const rp, 
		Real * const vm, Real * const vp,
		Real * const pm, Real * const pp,
		Real * const Gm, Real * const Gp, 
		Real * const PIm, Real * const PIp,
		Real * const outm, Real * const outp)
{
	enum { P = TempSOA::PITCH };

#define ID (ix + P*iy)

	for(int iy=0; iy<TempSOA::NY; ++iy)
		for(int ix=0; ix<NX; ix+=4)
		{
			const vector4double rminus = vec_ld(0L, rm + ID);
			const vector4double vminus = vec_ld(0L, vm + ID);
			const vector4double pminus = vec_ld(0L, pm + ID);
			const vector4double Gminus = vec_ld(0L, Gm + ID);
			const vector4double PIminus = vec_ld(0L, PIm + ID);

			const vector4double invGm = myreciprocal<preclevel>(Gminus);
			const vector4double cminus2 = vec_mul(vec_madd(invGm, vec_add(pminus, PIminus), pminus), myreciprocal<preclevel>(rminus));
			const vector4double cminus = mysqrt<preclevel>(cminus2);

			const vector4double rplus = vec_ld(0L, rp + ID);
			const vector4double vplus = vec_ld(0L, vp + ID);
			const vector4double pplus = vec_ld(0L, pp + ID);
			const vector4double Gplus = vec_ld(0L, Gp + ID);
			const vector4double PIplus = vec_ld(0L, PIp + ID);

			const vector4double invGp = myreciprocal<preclevel>(Gplus);
			const vector4double cplus2 = vec_mul(vec_madd(invGp, vec_add(pplus, PIplus), pplus), myreciprocal<preclevel>(rplus));
			const vector4double cplus = mysqrt<preclevel>(cplus2);

			const vector4double resultm = mymin(vec_sub(vminus, cminus), vec_sub(vplus, cplus)); 
			const vector4double resultp = mymax(vec_add(vminus, cminus), vec_add(vplus, cplus)); 

			vec_sta(resultm, 0L, outm + ID);
			vec_sta(resultp, 0L, outp + ID);
		}
#undef ID
}

void HLLESOA2D_QPX::rho(const TempSOA& rm, const TempSOA& rp,
		const TempSOA& vm, const TempSOA& vp,
		const TempSOA& am, const TempSOA& ap,
		TempSOA& out) const
{
	_qpx_hlle_rho<TempSOA::NX>(const_cast<Real*>(rm.ptr(0,0)), 
			const_cast<Real*>(rp.ptr(0,0)), 
			const_cast<Real*>(vm.ptr(0,0)), 
			const_cast<Real*>(vp.ptr(0,0)), 
			const_cast<Real*>(am.ptr(0,0)), 
			const_cast<Real*>(ap.ptr(0,0)), &out.ref(0,0));
}


void HLLESOA2D_QPX::vel(const TempSOA& rminus, const TempSOA& rplus,
		const TempSOA& vminus, const TempSOA& vplus,
		const TempSOA& vdminus, const TempSOA& vdplus,
		const TempSOA& aminus, const TempSOA& aplus,
		TempSOA& out) const
{
	_qpx_hlle_vel<TempSOA::NX>(const_cast<Real*>(rminus.ptr(0,0)), 
			const_cast<Real*>(rplus.ptr(0,0)), 
			const_cast<Real*>(vminus.ptr(0,0)), 
			const_cast<Real*>(vplus.ptr(0,0)),
			const_cast<Real*>(vdminus.ptr(0,0)), 
			const_cast<Real*>(vdplus.ptr(0,0)), 
			const_cast<Real*>(aminus.ptr(0,0)), 
			const_cast<Real*>(aplus.ptr(0,0)), &out.ref(0,0));
}


void HLLESOA2D_QPX::pvel(const TempSOA& rminus, const TempSOA& rplus,
		const TempSOA& vminus, const TempSOA& vplus,
		const TempSOA& pminus, const TempSOA& pplus,
		const TempSOA& aminus, const TempSOA& aplus,
		TempSOA& out) const
{
	_qpx_hlle_pvel<TempSOA::NX>(const_cast<Real*>(rminus.ptr(0,0)), 
			const_cast<Real*>(rplus.ptr(0,0)), 
			const_cast<Real*>(vminus.ptr(0,0)), 
			const_cast<Real*>(vplus.ptr(0,0)), 
			const_cast<Real*>(pminus.ptr(0,0)), 
			const_cast<Real*>(pplus.ptr(0,0)), 
			const_cast<Real*>(aminus.ptr(0,0)), 
			const_cast<Real*>(aplus.ptr(0,0)), &out.ref(0,0));
}


void HLLESOA2D_QPX::e(const TempSOA& rminus, const TempSOA& rplus,
		const TempSOA& vdminus, const TempSOA& vdplus,
		const TempSOA& v1minus, const TempSOA& v1plus,
		const TempSOA& v2minus, const TempSOA& v2plus,
		const TempSOA& pminus, const TempSOA& pplus,
		const TempSOA& Gminus, const TempSOA& Gplus, 
		const TempSOA& PIminus, const TempSOA& PIplus,
		const TempSOA& aminus, const TempSOA& aplus,
		TempSOA& out) const
{
	_qpx_hlle_e<TempSOA::NX>(const_cast<Real*>(rminus.ptr(0,0)), const_cast<Real*>(rplus.ptr(0,0)), 
			const_cast<Real*>(vdminus.ptr(0,0)), const_cast<Real*>(vdplus.ptr(0,0)), 
			const_cast<Real*>(v1minus.ptr(0,0)), const_cast<Real*>(v1plus.ptr(0,0)),
			const_cast<Real*>(v2minus.ptr(0,0)), const_cast<Real*>(v2plus.ptr(0,0)), 
			const_cast<Real*>(pminus.ptr(0,0)), const_cast<Real*>(pplus.ptr(0,0)), 
			const_cast<Real*>(Gminus.ptr(0,0)), const_cast<Real*>(Gplus.ptr(0,0)), 
			const_cast<Real*>(PIminus.ptr(0,0)), const_cast<Real*>(PIplus.ptr(0,0)),
			const_cast<Real*>(aminus.ptr(0,0)), const_cast<Real*>(aplus.ptr(0,0)), &out.ref(0,0));
}


void HLLESOA2D_QPX::char_vel(const TempSOA& rminus, const TempSOA& rplus, 
		const TempSOA& vminus, const TempSOA& vplus,
		const TempSOA& pminus, const TempSOA& pplus,
		const TempSOA& Gminus, const TempSOA& Gplus, 
		const TempSOA& Pminus, const TempSOA& Pplus,
		TempSOA& out_minus, TempSOA& out_plus) const
{
	_qpx_char_vel<TempSOA::NX>(const_cast<Real*>(rminus.ptr(0,0)), const_cast<Real*>(rplus.ptr(0,0)), 
			const_cast<Real*>(vminus.ptr(0,0)), const_cast<Real*>(vplus.ptr(0,0)), 
			const_cast<Real*>(pminus.ptr(0,0)), const_cast<Real*>(pplus.ptr(0,0)), 
			const_cast<Real*>(Gminus.ptr(0,0)), const_cast<Real*>(Gplus.ptr(0,0)), 
			const_cast<Real*>(Pminus.ptr(0,0)), const_cast<Real*>(Pplus.ptr(0,0)),
			&out_minus.ref(0,0), &out_plus.ref(0,0));
}

