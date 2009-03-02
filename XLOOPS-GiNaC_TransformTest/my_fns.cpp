/*******************************************************************************
**
**	xloop-ginacs Project
**	1Loop4Pt implementation
**
**
**	HCMUNS,
**	June 2008
**
**
**	Author(s): 	Son, D. H. 	(sondo01@gmail.com)
**			Khiem, Phan 	(phanhongkhiem@gmail.com)
**			Quan, Phan 	(anhquan.phanle@gmail.com)
********************************************************************************
**
**	Personal use functions
**
********************************************************************************
**
** Historial Log:
**	Date		Version	Author		Description
**	_________	_______	_________	________________________________
**	20090214	1.0	Quan Phan	Create this file
*******************************************************************************/

#include "ginac/ginac.h"
#include "my_fns.h"
#include "extrm.h"

using namespace GiNaC;

namespace xloops{
	/*20090211: Quan added this function to substitue the Rho and evaluation*/
	ex my_evalf(const ex &x){
		ex return_value;
	
	//	1.	Subs Rho's
		return_value = x.subs(lst(Rho == input_Rho, Rho1 == input_Rho1, Rho2 == input_Rho2));
	
	//	Last: evaluate
		return_value = return_value.evalf();
	
		return return_value;
	}
		
	ex my_is_zero_eval(const ex &x){
		ex factor = x;
		if(is_a<numeric>(factor)){
			if (abs(factor)<=1e-30)
				return 1.0;
			return 0.0;
		}
		return my_is_zero(x).hold();
	}
	REGISTER_FUNCTION(my_is_zero, eval_func(my_is_zero_eval));

	ex my_csgn_eval(const ex &x){
		ex factor = x;
		
		if (is_a<numeric>(factor)){
			if (my_is_zero(factor) == 1)
				return 0.0;
			return csgn(x);
		}
		return my_csgn(x).hold();
	}
	REGISTER_FUNCTION(my_csgn, eval_func(my_csgn_eval));

	ex my_step_eval(const ex &x){
		ex factor = x;
		
		if (is_a<numeric>(x)){
			if (real_part(x) >= 0.0)
				return 1.0;
			else
				return 0.0;
		}
		return my_step(x).hold();
	}
	REGISTER_FUNCTION(my_step, eval_func(my_step_eval));


	ex my_is_positive_eval(const ex &x){
		ex factor = x;
		
		if (is_a<numeric>(factor)) {
			if (factor.info(info_flags::positive))
				return 1.0;
			else
				return 0.0;
		}
		return my_step(x).hold();
	}
	REGISTER_FUNCTION(my_is_positive, eval_func(my_is_positive_eval));

	ex my_is_negative_eval(const ex &x){
		ex factor = x;
		
		if (is_a<numeric>(factor)) {
			if (factor.info(info_flags::negative))
				return 1.0;
			else
				return 0.0;
		}
		return my_step(x).hold();
	}
	REGISTER_FUNCTION(my_is_negative, eval_func(my_is_negative_eval));


	ex fn_eta_plus(const ex &sigma, const ex &z1, const ex &z2, const ex &P){
		ex eta_plus;
		eta_plus = my_step(imag_part(sigma * z1)) * my_step(imag_part(z2)) * my_step(P)
				-
				my_step(-imag_part(sigma*z1)) * my_step(-imag_part(z2)) * my_step(-P);
		eta_plus *= 2.0 * Pi * I;
		return eta_plus;
	}

	ex fn_eta_minus(const ex &sigma, const ex &z1, const ex &z2, const ex &P){
		ex eta_minus;
		eta_minus = my_step(-imag_part(sigma * z1)) * my_step(imag_part(z2)) * my_step(-P)
				-
				my_step(imag_part(sigma*z1)) * my_step(-imag_part(z2)) * my_step(P);
		eta_minus *= 2.0 * Pi * I;
		return eta_minus;
	}

	//	II.9	Misc. Functions
	ex fn_delta_eval (const ex & x){	// Delta Chroenacker
		if (is_a<numeric>(x))
			return my_is_zero(x);
		else
			return myfn_delta(x).hold();
	}

	REGISTER_FUNCTION(myfn_delta, eval_func(fn_delta_eval));

	ex myfn_IvZero_eval(const ex & x){
		if (my_is_zero(x) == 1.0)
			return 0.0;
		else if (my_is_zero(x) == 0.0)
			return 1.0;
		return myfn_IvZero(x).hold();
	}

	REGISTER_FUNCTION(myfn_IvZero, eval_func(myfn_IvZero_eval));

	ex myfn_Zero_eval(const ex & x){
		if (is_a<numeric>(x))
			return my_is_zero(x);
		else
			return myfn_Zero(x).hold();
	}

	REGISTER_FUNCTION(myfn_Zero, eval_func(myfn_Zero_eval));

	ex myfn_IvTheta_eval(const ex & x){
		if (is_a<numeric>(x)){
			if (real_part(x).info(info_flags::negative))
				return 1.0;
			else
				return 0.0;
		}
		else
			return myfn_IvTheta(x).hold();
	}

	REGISTER_FUNCTION(myfn_IvTheta, eval_func(myfn_IvTheta_eval));

	ex myfn_Sign0_eval(const ex & x){
		if (is_a<numeric>(x)){
			if (real_part(x).info(info_flags::negative))
				return -1.0;
			else
				return 1.0;
		}
		else
			return myfn_Sign0(x).hold();
	}

	REGISTER_FUNCTION(myfn_Sign0, eval_func(myfn_Sign0_eval));


	ex fn_eta_eval(const ex &x, const ex &y){
	 /** *****************************************************************************
		 **	eta(x,y) = 2pi.i.{
		 **			theta(-Im(x)).theta(-Im(y)).theta(Im(x.y))
		 **			- 
		 **			theta(Im(x)).theta(Im(y)).theta(-Im(x.y))
		 **			}
	  *********************************************************************************/
		if (is_a<numeric>(x) && is_a<numeric>(y)){
			ex eta;
			eta = 2.0*Pi*I*(
					my_step(-imag_part(x))*my_step(-imag_part(y))*my_step(imag_part(x*y))
					-
					my_step(imag_part(x))*my_step(imag_part(y))*my_step(-imag_part(x*y))
				       );
			return eta;
		}
		return fn_eta(x,y).hold();
	}
	REGISTER_FUNCTION(fn_eta, eval_func(fn_eta_eval));

	ex fn_R(const ex &x, const ex &y){
	 /** *****************************************************************************
		 **		ln(x) - ln(y)
		 **	R = ----------------------
		 **		  x - y
	  *********************************************************************************/
		if(x == y)
			throw std::runtime_error("R Function's Denominator = 0");
	
		return (log(x) - log(y)) / (x - y);
	}

	ex fn_GPositive(const ex &T1, const ex &T2){
	 /** *****************************************************************************
		 **	G+ = R(-T1, -T2)
	  *********************************************************************************/
		return fn_R(-T1, -T2);
	}

	ex fn_GNegative(const ex &T1, const ex &T2){
	 /** *****************************************************************************
		 **	G- = R(T1, T2)
	  *********************************************************************************/
		return fn_R(T1, T2);
	}

	ex fn_LogACG(const ex &a, const ex &b, const ex &x, const ex &y){
		ex LogACG_abxy;
		ex t, A, x0, y0, z0;
	
		if(a.info(info_flags::negative) || my_is_zero(a)==1.0)
			throw std::runtime_error("The a factor of LogACG isn't positive!");		
		/*step10calculation, eq 7*/
		t = b/a;
		A = sqrt(real_part(t) + imag_part(t)) + sqrt(real_part(x) + imag_part(x)) + sqrt(real_part(y) + imag_part(y));
		x0 = x/A; y0 = y/A; z0 = t/A;
	
		LogACG_abxy = (log(x0) - log(y0)) * log(a*A);
		LogACG_abxy -= pow(log(x0), 2)/2.0 + pow(log(y0), 2)/2.0 + Li2(1.0 - z0/y0) - Li2(1.0 - z0/x0)
				+
				log(y0)*( fn_eta(z0-y0, 1.0/(1.0 - y0)) - fn_eta(z0-y0, -1.0/y0))
				-
				log(x0)*( fn_eta(z0-x0, 1.0/(1.0 - x0)) - fn_eta(z0-y0, -1.0/x0))
				+
				log(1.0 - z0/y0) * fn_eta(z0, 1.0/y0)
				-
				log(1.0 - z0/x0) * fn_eta(z0, 1.0/x0)
				;
		LogACG_abxy /= A*(x0-y0);
		return LogACG_abxy;
	}

	ex fn_LogARG(const ex &a, const ex &b, const ex &x, const ex &y){
		ex LogARG_abxy;
		if(a.info(info_flags::positive) || my_is_zero(a)==1.0)
			throw std::runtime_error("The a factor of LogARG isn't negative!");
	
		LogARG_abxy = log(b) * (log(x) - log(y)) / (x - y) + fn_LogACG(a/b, 1.0, x, y);
	
		return LogARG_abxy;
	}

	ex fn_LogAG(const ex &a, const ex &b, const ex &x, const ex &y){
		ex LogAG_abxy;
		if(my_is_zero(a)==1.0)
			throw std::runtime_error("The a factor of LogARG is zero!");
	
		if(a.info(info_flags::positive))
			LogAG_abxy = fn_LogACG(a,b,x,y);
		if(a.info(info_flags::negative))
			LogAG_abxy = fn_LogARG(a,b,x,y);
		return LogAG_abxy;
	}

	ex p2q(const ex &p){
		ex p12, p23, p13, p1s, p2s, p3s, p4s, p12s, p23s;
		ex q, q10, q20, q21, q30, q31, q32;
		
		p1s = p.op(0); p2s = p.op(1); p3s = p.op(2); p4s = p.op(3);
		p12s = p.op(4); p23s = p.op(5);
			
		p12 = (p12s - p1s - p2s)/2.0;
		p23 = (p23s - p2s - p3s)/2.0;
		p13 = (p2s + p4s - p12s - p23s)/2.0;
		
		q10 = sqrt(p1s);
		q20 = p12/sqrt(p1s) + sqrt(p1s);
		q21 = sqrt( pow(p12/sqrt(p1s) + sqrt(p1s),2) - p12s );
		q30 = p13/sqrt(p1s) + p12/sqrt(p1s) + sqrt(p1s);
		q31 = (q30*q20 - p12s - p13 - p23) / q21;
		q32 = sqrt( pow(q30,2) - pow(q31,2) - p4s );
		
		q = lst(q10, q20, q21, q30, q31, q32);
		
		return q;
	}

}// Namespace xloops
