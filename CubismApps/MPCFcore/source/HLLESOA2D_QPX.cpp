#include "HLLESOA2D_QPX.h"
#include "common.h"

#ifndef _MICROFUSION_
#define _MICROFUSION_ 2
#endif

//=============================== KERNELS ===========================================

inline void _qpx_hlle_rho(Real * const rm, Real * const rp,
		Real * const vm, Real * const vp,
		Real * const am, Real * const ap,
		Real * const out)
{
	enum { NTOTAL = TempSOA::PITCH * TempSOA::NY};

	for(int ID = 0; ID < NTOTAL; ID += 4) 
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
}

inline void _qpx_hlle_vel(Real * const rm, Real * const rp,
		Real * const vm, Real * const vp,
		Real * const vdm, Real * const vdp,
		Real * const am, Real * const ap,
		Real * const out)
{
	enum { NTOTAL = TempSOA::PITCH * TempSOA::NY};

	for(int ID = 0; ID < NTOTAL; ID += 4) 
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
}


inline void _qpx_hlle_pvel(Real * const rm, Real * const rp,
		Real * const vm, Real * const vp,
		Real * const pm, Real * const pp,
		Real * const am, Real * const ap,
		Real * const out)
{
	enum { NTOTAL = TempSOA::PITCH * TempSOA::NY};

	for(int ID = 0; ID < NTOTAL; ID += 4) 
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
}

inline void _qpx_hlle_pvel_rho(Real * const rm, Real * const rp,
		Real * const vm, Real * const vp,
		Real * const pm, Real * const pp,
		Real * const am, Real * const ap,
		Real * const outrho, Real * const outpvel)
{
	enum { NTOTAL = TempSOA::PITCH * TempSOA::NY};

	for(int ID = 0; ID < NTOTAL; ID += 4)
	{
		const vector4double rminus = vec_lda(0L, rm + ID);
		const vector4double rpluss = vec_lda(0L, rp + ID);

		const vector4double myvminus = vec_lda(0L, vm + ID);
		const vector4double myvpluss = vec_lda(0L, vp + ID);

		const vector4double uminus = vec_mul(myvminus, rminus); 
		const vector4double upluss = vec_mul(myvpluss, rpluss);

		const vector4double aminus = vec_lda(0L, am + ID);
		const vector4double apluss = vec_lda(0L, ap + ID);

		const vector4double fminus = vec_madd(myvminus, uminus, vec_lda(0L, pm + ID));
		const vector4double fpluss = vec_madd(myvpluss, upluss, vec_lda(0L, pp + ID)); 

		const vector4double tmp1_rho = vec_mul(vec_mul(aminus, apluss), vec_sub(rpluss, rminus));
		const vector4double tmp1 = vec_mul(vec_mul(aminus, apluss), vec_sub(upluss, uminus));
		const vector4double tmp2_rho = vec_madd(apluss, uminus, vec_nmsub(aminus, upluss, tmp1_rho));
		const vector4double tmp2 = vec_madd(apluss, fminus, vec_nmsub(aminus, fpluss, tmp1));

		const vector4double fother_rho = vec_mul(tmp2_rho, myreciprocal<preclevel>(vec_sub(apluss, aminus)));
		const vector4double fother = vec_mul(tmp2, myreciprocal<preclevel>(vec_sub(apluss, aminus)));

		const vector4double result_rho = vec_sel(upluss, vec_sel(fother_rho, uminus, aminus), apluss);
		const vector4double result = vec_sel(fpluss, vec_sel(fother, fminus, aminus), apluss);

		vec_sta(result_rho, 0L, outrho + ID);
		vec_sta(result, 0L, outpvel + ID);
	}
}


inline void _qpx_hlle_e(Real * const rm, Real * const rp,
		Real * const vdm, Real * const vdp,
		Real * const v1m, Real * const v1p,
		Real * const v2m, Real * const v2p,
		Real * const pm, Real * const pp,
		Real * const Gm, Real * const Gp, 
		Real * const PIm, Real * const PIp,
		Real * const am, Real * const ap,
		Real * const out)
{

	const vector4double F_1_2 = vec_splats(0.5);

	enum { NTOTAL = TempSOA::PITCH * TempSOA::NY};

	for(int ID = 0; ID < NTOTAL; ID += 4) 
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
}

inline void _qpx_char_vel(Real * const rm, Real * const rp, 
		Real * const vm, Real * const vp,
		Real * const pm, Real * const pp,
		Real * const Gm, Real * const Gp, 
		Real * const PIm, Real * const PIp,
		Real * const outm, Real * const outp)
{
	enum { NTOTAL = TempSOA::PITCH * TempSOA::NY};

	for(int ID = 0; ID < NTOTAL; ID += 4) 
	{
		const vector4double rminus = vec_lda(0L, rm + ID);
		const vector4double vminus = vec_lda(0L, vm + ID);
		const vector4double pminus = vec_lda(0L, pm + ID);
		const vector4double Gminus = vec_lda(0L, Gm + ID);
		const vector4double PIminus = vec_lda(0L, PIm + ID);

		const vector4double invGm = myreciprocal<preclevel>(Gminus);
		const vector4double cminus2 = vec_mul(vec_madd(invGm, vec_add(pminus, PIminus), pminus), myreciprocal<preclevel>(rminus));
		const vector4double cminus = mysqrt<preclevel>(cminus2);

		const vector4double rplus = vec_lda(0L, rp + ID);
		const vector4double vplus = vec_lda(0L, vp + ID);
		const vector4double pplus = vec_lda(0L, pp + ID);
		const vector4double Gplus = vec_lda(0L, Gp + ID);
		const vector4double PIplus = vec_lda(0L, PIp + ID);

		const vector4double invGp = myreciprocal<preclevel>(Gplus);
		const vector4double cplus2 = vec_mul(vec_madd(invGp, vec_add(pplus, PIplus), pplus), myreciprocal<preclevel>(rplus));
		const vector4double cplus = mysqrt<preclevel>(cplus2);

		const vector4double resultm = mymin(vec_sub(vminus, cminus), vec_sub(vplus, cplus)); 
		const vector4double resultp = mymax(vec_add(vminus, cminus), vec_add(vplus, cplus)); 

		vec_sta(resultm, 0L, outm + ID);
		vec_sta(resultp, 0L, outp + ID);
	}
}

inline void _qpx_hllee_charvel(Real * const rm, Real * const rp,
		Real * const vdm, Real * const vdp,
		Real * const v1m, Real * const v1p,
		Real * const v2m, Real * const v2p,
		Real * const pm, Real * const pp,
		Real * const Gm, Real * const Gp,
		Real * const PIm, Real * const PIp,
		Real * const am, Real * const ap,
		Real * const oute, Real * const outcm, Real * const outcp)
{

	const vector4double F_1_2 = vec_splats(0.5);

	enum { NTOTAL = TempSOA::PITCH * TempSOA::NY};

	for(int ID = 0; ID < NTOTAL; ID += 4)
	{
		const vector4double rminus      = vec_lda(0L, rm + ID);
		const vector4double vdminus     = vec_lda(0L, vdm + ID);
		const vector4double v1minus     = vec_lda(0L, v1m + ID);
		const vector4double v2minus     = vec_lda(0L, v2m + ID);
		const vector4double pminus      = vec_lda(0L, pm + ID);
		const vector4double Gminus      = vec_lda(0L, Gm + ID);
		const vector4double PIminus     = vec_lda(0L, PIm + ID);

		const vector4double speedminus = vec_madd(vdminus, vdminus, vec_madd(v1minus, v1minus, vec_mul(v2minus, v2minus)));
		const vector4double eminus = vec_madd(pminus, Gminus, vec_madd(vec_mul(F_1_2, rminus), speedminus, PIminus));
		const vector4double invGm = myreciprocal<preclevel>(Gminus);
		const vector4double cminus2 = vec_mul(vec_madd(invGm, vec_add(pminus, PIminus), pminus), myreciprocal<preclevel>(rminus));
		const vector4double cminus = mysqrt<preclevel>(cminus2);

		const vector4double rplus       = vec_lda(0L, rp + ID);
		const vector4double vdplus      = vec_lda(0L, vdp + ID);
		const vector4double v1plus      = vec_lda(0L, v1p + ID);
		const vector4double v2plus      = vec_lda(0L, v2p + ID);
		const vector4double pplus       = vec_lda(0L, pp + ID);
		const vector4double Gplus       = vec_lda(0L, Gp + ID);
		const vector4double PIplus      = vec_lda(0L, PIp + ID);

		const vector4double speedplus = vec_madd(vdplus, vdplus, vec_madd(v1plus, v1plus, vec_mul(v2plus, v2plus)));
		const vector4double eplus = vec_madd(pplus, Gplus, vec_madd(vec_mul(F_1_2, rplus), speedplus, PIplus));
		const vector4double invGp = myreciprocal<preclevel>(Gplus);
		const vector4double cplus2 = vec_mul(vec_madd(invGp, vec_add(pplus, PIplus), pplus), myreciprocal<preclevel>(rplus));
		const vector4double cplus = mysqrt<preclevel>(cplus2);

		const vector4double aminus      = mymin(vec_sub(vdminus, cminus), vec_sub(vdplus, cplus));
		const vector4double apluss      = mymax(vec_add(vdminus, cminus), vec_add(vdplus, cplus));
		const vector4double fminus = vec_mul(vdminus, vec_add(pminus, eminus));
		const vector4double fpluss = vec_mul(vdplus, vec_add(pplus, eplus));

		const vector4double tmp1 = vec_mul(vec_mul(aminus, apluss), vec_sub(eplus, eminus));
		const vector4double tmp2 = vec_madd(apluss, fminus, vec_nmsub(aminus, fpluss, tmp1));
		const vector4double fother = vec_mul(tmp2, myreciprocal<preclevel>(vec_sub(apluss, aminus)));
		const vector4double result = vec_sel(fpluss, vec_sel(fother, fminus, aminus), apluss);

		vec_sta(result, 0L, oute + ID);
		vec_sta(aminus, 0L, outcm + ID);
		vec_sta(apluss, 0L, outcp + ID);
	}
}

inline void _qpx_hlle_all(Real * const rm, Real * const rp,
						  Real * const vdm, Real * const vdp,
						  Real * const v1m, Real * const v1p,
						  Real * const v2m, Real * const v2p,
						  Real * const pm, Real * const pp,
						  Real * const Gm, Real * const Gp, 
						  Real * const PIm, Real * const PIp,
						  Real * const outam, Real * const outap, Real * const outrho, 
						  Real * const outvd, Real * const outv1, Real * const outv2,
						  Real * const oute, Real * const outG, Real * const outP)
{
	enum { NTOTAL = TempSOA::PITCH * TempSOA::NY};
	
	const vector4double F_1_2 = vec_splats(0.5);
	
	for(int ID = 0; ID < NTOTAL; ID += 4) 
	{
		const vector4double rminus = vec_lda(0L, rm + ID);
		const vector4double vdminus = vec_lda(0L, vdm + ID);
		const vector4double v1minus	= vec_lda(0L, v1m + ID);
		const vector4double v2minus	= vec_lda(0L, v2m + ID);
		const vector4double pminus = vec_lda(0L, pm + ID);
		const vector4double Gminus = vec_lda(0L, Gm + ID);
		const vector4double PIminus = vec_lda(0L, PIm + ID);
		
		const vector4double uminus = vec_mul(vdminus, rminus);
		const vector4double uminus_v1 = vec_mul(v1minus, rminus);
		const vector4double uminus_v2 = vec_mul(v2minus, rminus);
		const vector4double speedminus = vec_madd(vdminus, vdminus, vec_madd(v1minus, v1minus, vec_mul(v2minus, v2minus)));
		const vector4double eminus = vec_madd(pminus, Gminus, vec_madd(vec_mul(F_1_2, rminus), speedminus, PIminus));
		const vector4double invGm = myreciprocal<preclevel>(Gminus);
		const vector4double cminus2 = vec_mul(vec_madd(invGm, vec_add(pminus, PIminus), pminus), myreciprocal<preclevel>(rminus));
		const vector4double cminus = mysqrt<preclevel>(cminus2);
		
		const vector4double rplus = vec_lda(0L, rp + ID);
		const vector4double vdplus = vec_lda(0L, vdp + ID);
		const vector4double v1plus = vec_lda(0L, v1p + ID);
		const vector4double v2plus = vec_lda(0L, v2p + ID);
		const vector4double pplus = vec_lda(0L, pp + ID);
		const vector4double Gplus = vec_lda(0L, Gp + ID);
		const vector4double PIplus = vec_lda(0L, PIp + ID);
		
		const vector4double uplus = vec_mul(vdplus, rplus);
		const vector4double uplus_v1 = vec_mul(v1plus, rplus);
		const vector4double uplus_v2 = vec_mul(v2plus, rplus);
		const vector4double speedplus = vec_madd(vdplus, vdplus, vec_madd(v1plus, v1plus, vec_mul(v2plus, v2plus)));
		const vector4double eplus = vec_madd(pplus, Gplus, vec_madd(vec_mul(F_1_2, rplus), speedplus, PIplus));
		const vector4double invGp = myreciprocal<preclevel>(Gplus);
		const vector4double cplus2 = vec_mul(vec_madd(invGp, vec_add(pplus, PIplus), pplus), myreciprocal<preclevel>(rplus));
		const vector4double cplus = mysqrt<preclevel>(cplus2);
		
		const vector4double aminus = mymin(vec_sub(vdminus, cminus), vec_sub(vdplus, cplus)); 
		const vector4double aplus = mymax(vec_add(vdminus, cminus), vec_add(vdplus, cplus)); 
		vec_sta(aminus, 0L, outam + ID);
		vec_sta(aplus, 0L, outap + ID);
		
		const vector4double fminus_rho = vec_mul(vdminus, rminus);
		const vector4double fminus_vd = vec_madd(vdminus, uminus, pminus);
		const vector4double fminus_v1 = vec_mul(vdminus, uminus_v1);
		const vector4double fminus_v2 = vec_mul(vdminus, uminus_v2);
		const vector4double fminus_e = vec_mul(vdminus, vec_add(pminus, eminus));
		const vector4double fminus_G = vec_mul(vdminus, Gminus);
		const vector4double fminus_P = vec_mul(vdminus, PIminus);
		
		const vector4double fplus_rho = vec_mul(vdplus, rplus);
		const vector4double fplus_vd = vec_madd(vdplus, uplus, pplus);
		const vector4double fplus_v1 = vec_mul(vdplus, uplus_v1);
		const vector4double fplus_v2 = vec_mul(vdplus, uplus_v2);
		const vector4double fplus_e = vec_mul(vdplus, vec_add(pplus, eplus));
		const vector4double fplus_G = vec_mul(vdplus, Gplus);
		const vector4double fplus_P = vec_mul(vdplus, PIplus);
		
		const vector4double amul = vec_mul(aminus, aplus);
		
		const vector4double tmp1_rho = vec_mul(amul, vec_sub(rplus, rminus));
		const vector4double tmp1_vd = vec_mul(amul, vec_sub(uplus, uminus));
		const vector4double tmp1_v1 = vec_mul(amul, vec_sub(uplus_v1, uminus_v1));
		const vector4double tmp1_v2 = vec_mul(amul, vec_sub(uplus_v2, uminus_v2));
		const vector4double tmp1_e = vec_mul(amul, vec_sub(eplus, eminus));
		const vector4double tmp1_G = vec_mul(amul, vec_sub(Gplus, Gminus));
		const vector4double tmp1_P = vec_mul(amul, vec_sub(PIplus, PIminus));
		
		const vector4double tmp2_rho = vec_madd(aplus, fminus_rho, vec_nmsub(aminus, fplus_rho, tmp1_rho));
		const vector4double tmp2_vd = vec_madd(aplus, fminus_vd, vec_nmsub(aminus, fplus_vd, tmp1_vd));
		const vector4double tmp2_v1 = vec_madd(aplus, fminus_v1, vec_nmsub(aminus, fplus_v1, tmp1_v1));
		const vector4double tmp2_v2 = vec_madd(aplus, fminus_v2, vec_nmsub(aminus, fplus_v2, tmp1_v2));
		const vector4double tmp2_e = vec_madd(aplus, fminus_e, vec_nmsub(aminus, fplus_e, tmp1_e));
		const vector4double tmp2_G = vec_madd(aplus, fminus_G, vec_nmsub(aminus, fplus_G, tmp1_G));
		const vector4double tmp2_P = vec_madd(aplus, fminus_P, vec_nmsub(aminus, fplus_P, tmp1_P));
		
		const vector4double inv_adiff = myreciprocal<preclevel>(vec_sub(aplus, aminus));
		
		const vector4double fother_rho = vec_mul(tmp2_rho, inv_adiff);
		const vector4double fother_vd = vec_mul(tmp2_vd, inv_adiff);
		const vector4double fother_v1 = vec_mul(tmp2_v1, inv_adiff);
		const vector4double fother_v2 = vec_mul(tmp2_v2, inv_adiff);
		const vector4double fother_e = vec_mul(tmp2_e, inv_adiff);
		const vector4double fother_G = vec_mul(tmp2_G, inv_adiff);
		const vector4double fother_P = vec_mul(tmp2_P, inv_adiff);
		
		const vector4double result_rho = vec_sel(fplus_rho, vec_sel(fother_rho, fminus_rho, aminus), aplus);
		const vector4double result_vd = vec_sel(fplus_vd, vec_sel(fother_vd, fminus_vd, aminus), aplus);
		const vector4double result_v1 = vec_sel(fplus_v1, vec_sel(fother_v1, fminus_v1, aminus), aplus);
		const vector4double result_v2 = vec_sel(fplus_v2, vec_sel(fother_v2, fminus_v2, aminus), aplus);
		const vector4double result_e = vec_sel(fplus_e, vec_sel(fother_e, fminus_e, aminus), aplus);
		const vector4double result_G = vec_sel(fplus_G, vec_sel(fother_G, fminus_G, aminus), aplus);
		const vector4double result_P = vec_sel(fplus_P, vec_sel(fother_P, fminus_P, aminus), aplus);
		
		vec_sta(result_rho, 0L, outrho + ID);
		vec_sta(result_vd, 0L, outvd + ID);
		vec_sta(result_v1, 0L, outv1 + ID);
		vec_sta(result_v2, 0L, outv2 + ID);
		vec_sta(result_e, 0L, oute + ID);
		vec_sta(result_G, 0L, outG + ID);
		vec_sta(result_P, 0L, outP + ID);
	}	
}

//=============================== DRIVERS ===========================================


void HLLESOA2D_QPX::rho(const TempSOA& rm, const TempSOA& rp,
		const TempSOA& vm, const TempSOA& vp,
		const TempSOA& am, const TempSOA& ap,
		TempSOA& out) const
{
	_qpx_hlle_rho(const_cast<Real*>(rm.ptr(0,0)), 
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
	_qpx_hlle_vel(const_cast<Real*>(rminus.ptr(0,0)), 
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
	_qpx_hlle_pvel(const_cast<Real*>(rminus.ptr(0,0)), 
			const_cast<Real*>(rplus.ptr(0,0)), 
			const_cast<Real*>(vminus.ptr(0,0)), 
			const_cast<Real*>(vplus.ptr(0,0)), 
			const_cast<Real*>(pminus.ptr(0,0)), 
			const_cast<Real*>(pplus.ptr(0,0)), 
			const_cast<Real*>(aminus.ptr(0,0)), 
			const_cast<Real*>(aplus.ptr(0,0)), &out.ref(0,0));
}


void HLLESOA2D_QPX::pvel_rho(const TempSOA& rminus, const TempSOA& rplus,
		const TempSOA& vminus, const TempSOA& vplus,
		const TempSOA& pminus, const TempSOA& pplus,
		const TempSOA& aminus, const TempSOA& aplus,
		TempSOA& outrho, TempSOA& outpvel ) const
{
	_qpx_hlle_pvel_rho(const_cast<Real*>(rminus.ptr(0,0)),
			const_cast<Real*>(rplus.ptr(0,0)),
			const_cast<Real*>(vminus.ptr(0,0)),
			const_cast<Real*>(vplus.ptr(0,0)),
			const_cast<Real*>(pminus.ptr(0,0)),
			const_cast<Real*>(pplus.ptr(0,0)),
			const_cast<Real*>(aminus.ptr(0,0)),
			const_cast<Real*>(aplus.ptr(0,0)), &outrho.ref(0,0), &outpvel.ref(0,0));
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
	_qpx_hlle_e(const_cast<Real*>(rminus.ptr(0,0)), const_cast<Real*>(rplus.ptr(0,0)), 
			const_cast<Real*>(vdminus.ptr(0,0)), const_cast<Real*>(vdplus.ptr(0,0)), 
			const_cast<Real*>(v1minus.ptr(0,0)), const_cast<Real*>(v1plus.ptr(0,0)),
			const_cast<Real*>(v2minus.ptr(0,0)), const_cast<Real*>(v2plus.ptr(0,0)), 
			const_cast<Real*>(pminus.ptr(0,0)), const_cast<Real*>(pplus.ptr(0,0)), 
			const_cast<Real*>(Gminus.ptr(0,0)), const_cast<Real*>(Gplus.ptr(0,0)), 
			const_cast<Real*>(PIminus.ptr(0,0)), const_cast<Real*>(PIplus.ptr(0,0)),
			const_cast<Real*>(aminus.ptr(0,0)), const_cast<Real*>(aplus.ptr(0,0)), &out.ref(0,0));
}

void HLLESOA2D_QPX::e_charvel(const TempSOA& rminus, const TempSOA& rplus,
		const TempSOA& vdminus, const TempSOA& vdplus,
		const TempSOA& v1minus, const TempSOA& v1plus,
		const TempSOA& v2minus, const TempSOA& v2plus,
		const TempSOA& pminus, const TempSOA& pplus,
		const TempSOA& Gminus, const TempSOA& Gplus,
		const TempSOA& PIminus, const TempSOA& PIplus,
		const TempSOA& aminus, const TempSOA& aplus,
		TempSOA& oute, TempSOA& outcm, TempSOA& outcp) const
{
	_qpx_hllee_charvel(const_cast<Real*>(rminus.ptr(0,0)), const_cast<Real*>(rplus.ptr(0,0)),
			const_cast<Real*>(vdminus.ptr(0,0)), const_cast<Real*>(vdplus.ptr(0,0)),
			const_cast<Real*>(v1minus.ptr(0,0)), const_cast<Real*>(v1plus.ptr(0,0)),
			const_cast<Real*>(v2minus.ptr(0,0)), const_cast<Real*>(v2plus.ptr(0,0)),
			const_cast<Real*>(pminus.ptr(0,0)), const_cast<Real*>(pplus.ptr(0,0)),
			const_cast<Real*>(Gminus.ptr(0,0)), const_cast<Real*>(Gplus.ptr(0,0)),
			const_cast<Real*>(PIminus.ptr(0,0)), const_cast<Real*>(PIplus.ptr(0,0)),
			const_cast<Real*>(aminus.ptr(0,0)), const_cast<Real*>(aplus.ptr(0,0)),
			&oute.ref(0,0), &outcm.ref(0,0), &outcp.ref(0,0));
}


void HLLESOA2D_QPX::char_vel(const TempSOA& rminus, const TempSOA& rplus, 
		const TempSOA& vminus, const TempSOA& vplus,
		const TempSOA& pminus, const TempSOA& pplus,
		const TempSOA& Gminus, const TempSOA& Gplus, 
		const TempSOA& Pminus, const TempSOA& Pplus,
		TempSOA& out_minus, TempSOA& out_plus) const
{
	_qpx_char_vel(const_cast<Real*>(rminus.ptr(0,0)), const_cast<Real*>(rplus.ptr(0,0)), 
			const_cast<Real*>(vminus.ptr(0,0)), const_cast<Real*>(vplus.ptr(0,0)), 
			const_cast<Real*>(pminus.ptr(0,0)), const_cast<Real*>(pplus.ptr(0,0)), 
			const_cast<Real*>(Gminus.ptr(0,0)), const_cast<Real*>(Gplus.ptr(0,0)), 
			const_cast<Real*>(Pminus.ptr(0,0)), const_cast<Real*>(Pplus.ptr(0,0)),
			&out_minus.ref(0,0), &out_plus.ref(0,0));
}


void HLLESOA2D_QPX::all(const TempSOA& rminus, const TempSOA& rplus,
		const TempSOA& vdminus, const TempSOA& vdplus,
		const TempSOA& v1minus, const TempSOA& v1plus,
		const TempSOA& v2minus, const TempSOA& v2plus,
		const TempSOA& pminus, const TempSOA& pplus,
		const TempSOA& Gminus, const TempSOA& Gplus,
		const TempSOA& PIminus, const TempSOA& PIplus,
		TempSOA& outam, TempSOA& outap, TempSOA& outrho, 
		TempSOA& outvd, TempSOA& outv1, TempSOA& outv2,
		TempSOA& oute,TempSOA& outG, TempSOA& outP) const
{
#if _MICROFUSION_ == 0
	
	char_vel(rminus, rplus, vdminus, vdplus, pminus, pplus, Gminus, Gplus, PIminus, PIplus, outam, outap);
	rho(rminus, rplus, vdminus, vdplus, outam, outap, outrho);
	pvel(rminus, rplus, vdminus, vdplus, pminus, pplus, outam, outap, outvd);
	vel(rminus, rplus, v1minus, v1plus, vdminus,  vdplus, outam, outap, outv1);
	vel(rminus, rplus, v2minus, v2plus, vdminus,  vdplus, outam, outap, outv2);
	e(rminus, rplus, vdminus, vdplus, v1minus, v1plus, v2minus, v2plus, pminus, pplus, Gminus, Gplus, PIminus, PIplus, outam, outap, oute);
	rho(Gminus, Gplus, vdminus, vdplus, outam, outap, outG);
	rho(PIminus, PIplus, vdminus, vdplus, outam, outap, outP);	
	
#elif _MICROFUSION_ == 1
	
	e_charvel(rminus, rplus, vdminus, vdplus, v1minus, v1plus, v2minus, v2plus, pminus, pplus, Gminus, Gplus, PIminus, PIplus, outam, outap, oute, outam, outap);
	pvel_rho(rminus, rplus, vdminus, vdplus, pminus, pplus, outam, outap, outrho, outvd);
	vel(rminus, rplus, v1minus, v1plus, vdminus,  vdplus, outam, outap, outv1);
	vel(rminus, rplus, v2minus, v2plus, vdminus,  vdplus, outam, outap, outv2);
	rho(Gminus, Gplus, vdminus, vdplus, outam, outap, outG);
	rho(PIminus, PIplus, vdminus, vdplus, outam, outap, outP);	
	
#else
	
	_qpx_hlle_all(const_cast<Real*>(rminus.ptr(0,0)), const_cast<Real*>(rplus.ptr(0,0)),
			const_cast<Real*>(vdminus.ptr(0,0)), const_cast<Real*>(vdplus.ptr(0,0)),
			const_cast<Real*>(v1minus.ptr(0,0)), const_cast<Real*>(v1plus.ptr(0,0)),
			const_cast<Real*>(v2minus.ptr(0,0)), const_cast<Real*>(v2plus.ptr(0,0)),
			const_cast<Real*>(pminus.ptr(0,0)), const_cast<Real*>(pplus.ptr(0,0)),
			const_cast<Real*>(Gminus.ptr(0,0)), const_cast<Real*>(Gplus.ptr(0,0)),
			const_cast<Real*>(PIminus.ptr(0,0)), const_cast<Real*>(PIplus.ptr(0,0)),
			&outam.ref(0,0), &outap.ref(0,0),  &outrho.ref(0,0), 
			&outvd.ref(0,0), &outv1.ref(0,0),  &outv2.ref(0,0),
			&oute.ref(0,0), &outG.ref(0,0),  &outP.ref(0,0));
#endif
}
