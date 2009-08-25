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
**	Level 2 functions
**
********************************************************************************
**
** Historial Log:
**	Date		Version	Author		Description
**	_________	_______	_________	________________________________
**	20090413	1.1	Quan Phan	Change the way to calculate.
								now can get the value directly.
**	20090214	1.0	Quan Phan	Create this file
*******************************************************************************/

#include "ginac/ginac.h"
#include "extrm.h"
#include "my_fns.h"
#include "trmchk.h"
#include "lev1.h"
#include "lev2.h"

using namespace GiNaC;
namespace xloops{

	ex fn_AC (int l, int k){
		// init
	  	ex AC_lk;
		ex a_lk = fn_a(l, k), c_lk = fn_c(l, k);
	
		// calc.
		AC_lk = a_lk + c_lk;
		return AC_lk;
	}

	ex fn_alpha (int l, int k){
		//init
		ex alpha_lk;
		ex AC_lk = fn_AC(l, k), b_lk = fn_a(l, k);
		
		check0denom(AC_lk, "alpha", l, k);
		
		// calc.
		alpha_lk = b_lk / AC_lk;
		
		return alpha_lk;
	}

	ex fn_A (int m, int l, int k){
		// init
		ex A_mlk;
		ex AC_mk = fn_AC(m, k), AC_lk = fn_AC(l, k), a_mk = fn_a(m, k), a_lk = fn_a(l, k);

		check0denom(AC_lk, "A", m, l, k);
		
		// calc.
		A_mlk = a_mk - a_lk*AC_mk/AC_lk;
		
		return A_mlk;
	}

	ex fn_B (int m, int l, int k){
		
		// init
		ex B_mlk;
		ex AC_mk = fn_AC(m, k), AC_lk = fn_AC(l, k), b_mk = fn_b(m, k), b_lk = fn_b(l, k);

		check0denom(AC_lk, "B", m, l, k);
		
		//calc.
		B_mlk = b_mk - b_lk*AC_mk/AC_lk;
		return B_mlk;
	}

	ex fn_C (int m, int l, int k){
		ex C_mlk;
		ex AC_mk = fn_AC(m, k), AC_lk = fn_AC(l, k), d_mk = fn_d(m, k), d_lk = fn_d(l, k);

		check0denom(AC_lk, "C", m, l, k);
	
		C_mlk = d_mk - d_lk*AC_mk/AC_lk;
		return C_mlk;
	}
	ex fn_C_im (int m, int l, int k){
		ex C_mlk_im;
	 	C_mlk_im = imag_part(fn_C(m, l, k));
		return C_mlk_im;
	}
	ex fn_C_re (int m, int l, int k){
		ex C_mlk_re;
	 	C_mlk_re = real_part(fn_C(m, l, k));
		return C_mlk_re;
	}
	ex fn_C_conj (int m, int l, int k){
		ex C_mlk_conj;
	 	C_mlk_conj = conjugate(fn_C(m, l, k));
		return C_mlk_conj;
	}

	ex fn_D (int m, int l, int k){
				// init
		ex a_lk = fn_a(l, k), AC_lk = fn_AC(l, k), alpha_lk = fn_alpha(l, k);

		ex D_mlk;
		
		// calc.
		D_mlk = 1.0 - 2.0*a_lk/AC_lk + pow(alpha_lk ,2);
		return D_mlk;
	}

	ex myfn_f_eval(const ex &factor){
		if(is_a<numeric>(factor)){
			if (my_is_zero(factor) == 1.0)
				return 1.0;
			if (factor.info(info_flags::negative))
				return 0.0;
			return 2.0;
		}
		else{
			return myfn_f(factor).hold();
		}
	}
	REGISTER_FUNCTION(myfn_f, eval_func(myfn_f_eval));

	ex fn_f (int l, int k){
		ex AC_lk = fn_AC(l, k), im_d_lk = fn_d_im(l, k);

		check0denom(AC_lk, "f", l, k);
		
	
		ex factor = - im_d_lk/AC_lk;
		return myfn_f(factor);
	}

	ex myfn_fminus_eval(const ex &factor){
		if(is_a<numeric>(factor)){
			if (my_is_zero(factor) == 1.0)
				return 1.0;
			if (factor.info(info_flags::negative))
				return 2.0;
			return 0.0;
		}
		else{
			return myfn_fminus(factor).hold();
		}
	}

	REGISTER_FUNCTION(myfn_fminus, eval_func(myfn_fminus_eval));
	ex fn_fminus (int l, int k){
		ex AC_lk = fn_AC(l, k), im_d_lk = fn_d_im(l, k);
		check0denom(AC_lk, "fminus", l, k);
		
		ex factor = - im_d_lk/AC_lk;
		return myfn_fminus(factor);
	}
}// Namespace xloops

