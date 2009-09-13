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
//		ex a_lk = fn_a(l, k), c_lk = fn_c(l, k);
		const ex a_lk = mat_a[l][k], c_lk = mat_c[l][k];
	
		// calc.
		AC_lk = a_lk + c_lk;
		return AC_lk;
	}

	ex fn_alpha (int l, int k){
		//init
		ex alpha_lk;
//		ex AC_lk = fn_AC(l, k), b_lk = fn_b(l, k);
		const ex AC_lk = mat_AC[l][k], b_lk = mat_b[l][k];
		
		check0denom(AC_lk, "alpha", l, k);
		
		// calc.
		alpha_lk = b_lk / AC_lk;
		
		return alpha_lk;
	}

	ex fn_A (int m, int l, int k){
		// init
		ex A_mlk;
//		ex AC_mk = fn_AC(m, k), AC_lk = fn_AC(l, k), a_mk = fn_a(m, k), a_lk = fn_a(l, k);
		const ex AC_mk = mat_AC[m][k], AC_lk = mat_AC[l][k], a_mk = mat_a[m][k], a_lk = mat_a[l][k];

		check0denom(AC_lk, "A", m, l, k);
		
		// calc.
		A_mlk = a_mk - a_lk*AC_mk/AC_lk;
		
		return A_mlk;
	}

	ex fn_B (int m, int l, int k){
		
		// init
		ex B_mlk;
//		ex AC_mk = fn_AC(m, k), AC_lk = fn_AC(l, k), b_mk = fn_b(m, k), b_lk = fn_b(l, k);
		const ex AC_mk = mat_AC[m][k], AC_lk = mat_AC[l][k], b_mk = mat_b[m][k], b_lk = mat_b[l][k];

		check0denom(AC_lk, "B", m, l, k);
		
		//calc.
		B_mlk = b_mk - b_lk*AC_mk/AC_lk;
		return B_mlk;
	}

	ex fn_C (int m, int l, int k){
		ex C_mlk;
//		ex AC_mk = fn_AC(m, k), AC_lk = fn_AC(l, k), d_mk = fn_d(m, k), d_lk = fn_d(l, k);
		const ex AC_mk = mat_AC[m][k], AC_lk = mat_AC[l][k], d_mk = mat_d[m][k], d_lk = mat_d[l][k];

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
//		ex a_lk = fn_a(l, k), AC_lk = fn_AC(l, k), alpha_lk = fn_alpha(l, k);
		const ex a_lk = mat_a[l][k], AC_lk = mat_AC[l][k], alpha_lk = mat_alpha[l][k];

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
//		ex AC_lk = fn_AC(l, k), im_d_lk = fn_d_im(l, k);
		const ex AC_lk = mat_AC[l][k], im_d_lk = mat_d_im[l][k];

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
//		ex AC_lk = fn_AC(l, k), im_d_lk = fn_d_im(l, k);
		const ex AC_lk = mat_AC[l][k], im_d_lk = mat_d_im[l][k];
		check0denom(AC_lk, "fminus", l, k);
		
		ex factor = - im_d_lk/AC_lk;
		return myfn_fminus(factor);
	}

	void lev2Calc(){
		// calculate previous level
		lev1Calc();

		int m, l, k;
		for (k = 0; k<4; k++) for (l=0; l<4; l++) if(l!=k){
			mat_AC[l][k] = fn_AC(l, k);
			mat_f[l][k] = fn_f(l, k);
			mat_fminus[l][k] = fn_fminus(l, k);
			mat_alpha[l][k] = fn_alpha(l, k);
		}
		for(k=0; k<4; k++) for(l=0; l<4; l++) for(m=0; m<4; m++) if(m!=l && m!=k && l!=k){
			mat_A[m][l][k] = fn_A(m, l, k);
			mat_B[m][l][k] = fn_B(m, l, k);
			mat_C[m][l][k] = fn_C(m, l, k);
			mat_C_im[m][l][k] = fn_C_im(m, l, k);
			mat_C_re[m][l][k] = fn_C_re(m, l, k);
			mat_C_conj[m][l][k] = fn_C_conj(m, l, k);
			mat_D[m][l][k] = fn_D(m, l, k);
		}
	}

	void lev2Calc(int debug){
		if (debug == 0){
			lev2Calc();
			return;
		}
		// calculate previous level
		lev1Calc(1);

		printf("Level 2 calculations...\n");

		int m, l, k;
		for (k = 0; k<4; k++) for (l=0; l<4; l++) if(l!=k){
			printf("AC[%d,%d]\n", l, k);
			mat_AC[l][k] = fn_AC(l, k);
			printf("f[%d,%d]\n", l, k);
			mat_f[l][k] = fn_f(l, k);
			printf("fminus[%d,%d]\n", l, k);
			mat_fminus[l][k] = fn_fminus(l, k);
			printf("alpha[%d,%d]\n", l, k);
			mat_alpha[l][k] = fn_alpha(l, k);
		}
		for(k=0; k<4; k++) for(l=0; l<4; l++) for(m=0; m<4; m++) if(m!=l && m!=k && l!=k){
			printf("A[%d,%d,%d]\n", m, l, k);
			mat_A[m][l][k] = fn_A(m, l, k);
			printf("B[%d,%d,%d]\n", m, l, k);
			mat_B[m][l][k] = fn_B(m, l, k);
			printf("C[%d,%d,%d]\n", m, l, k);
			mat_C[m][l][k] = fn_C(m, l, k);
			printf("C_im[%d,%d,%d]\n", m, l, k);
			mat_C_im[m][l][k] = fn_C_im(m, l, k);
			printf("C_re[%d,%d,%d]\n", m, l, k);
			mat_C_re[m][l][k] = fn_C_re(m, l, k);
			printf("C_conj[%d,%d,%d]\n", m, l, k);
			mat_C_conj[m][l][k] = fn_C_conj(m, l, k);
			printf("D[%d,%d,%d]\n", m, l, k);
			mat_D[m][l][k] = fn_D(m, l, k);
		}
	}

}// Namespace xloops

