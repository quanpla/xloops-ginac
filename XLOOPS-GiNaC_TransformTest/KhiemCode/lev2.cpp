/*******************************************************************************
**
**      xloop-ginacs Project
**      1Loop4Pt implementation
**
**
**      HCMUNS,
**      June 2008
**
**
**      Author(s):      Son, D. H.      (sondo01@gmail.com)
**      Khiem, Phan     (phanhongkhiem@gmail.com)
**      Quan, Phan      (anhquan.phanle@gmail.com)
********************************************************************************
**
**      Level 2 = 3 dimensions D0 integral
**
********************************************************************************
**
** Historial Log:
**      Date    	Version Author  Description
**      _________       _______ _________       ________________________________
**      20090221	1.0     Quan Phan       Create this file
*******************************************************************************/

#include "ginac/ginac.h"
#include "extrm.h" // for extern definitions of terms
#include "trmchk.h" // for checking the value of each term, for example, is a denominator zero?
#include "my_fns.h" // misc. function
#include "lev2.h"

using namespace GiNaC;
namespace xloops{

//      Level 2 Variable Functions
	ex fn_AC (int l, int k){
 /**
		 **     AC_lk = a_lk + c_lk
  **/
		ex AC_lk;
		AC_lk = mat_a[l][k] + mat_c[l][k];
		return AC_lk;
	}

	ex fn_alpha (int l, int k){
		ex alpha_lk;
 /**
		 **     b_lk
		 **     alpha_lk = -------------
		 **a_lk  +  c_lk
		 ** If denominator = 0, we set alpha_lk = 0 by default
  **/
		ex denom = mat_a[l][k] + mat_c[l][k];
       
		if(my_is_zero(denom) == 1.0){
			alpha_lk = 0.0;
		}
		else{
			alpha_lk = mat_b[l][k] / denom;
		}
		return alpha_lk;
	}
	ex fn_A (int m, int l, int k){
		ex A_mlk;
 /**
		 **  AC_mk
		 **     A_mlk = a_mk - a_lk -------
		 **  AC_lk
  **/
		ex denom = mat_AC[l][k];
		check0denom(denom, "A", m, l, k);

		A_mlk = mat_a[m][k] - mat_a[l][k]*mat_AC[m][k]/mat_AC[l][k];
		return A_mlk;
	}

	ex fn_B (int m, int l, int k){
		ex B_mlk;
 /**
		 **  AC_mk
		 **     B_mlk = b_mk - b_lk -------
		 **  AC_lk
  **/
		ex denom = mat_AC[l][k];
		check0denom(denom, "B", m, l, k);

		B_mlk = mat_b[m][k] - mat_b[l][k]*mat_AC[m][k]/mat_AC[l][k];
		return B_mlk;
	}

	ex fn_C (int m, int l, int k){
		ex C_mlk;
 /**
		 **  AC_mk
		 **     C_mlk = d_mk - d_lk -------
		 **  AC_lk
  **/
		ex denom = mat_AC[l][k];
		check0denom(denom, "C", m, l, k);

		C_mlk = mat_d[m][k] - mat_d[l][k]*mat_AC[m][k]/mat_AC[l][k];
		return C_mlk;
	}
	ex fn_C_im (int m, int l, int k){
		ex C_mlk_im;
		C_mlk_im = imag_part(mat_C[m][l][k]);
		return C_mlk_im;
	}
	ex fn_C_re (int m, int l, int k){
		ex C_mlk_re;
		C_mlk_re = real_part(mat_C[m][l][k]);
		return C_mlk_re;
	}
	ex fn_C_conj (int m, int l, int k){
		ex C_mlk_conj;
		C_mlk_conj = mat_C[m][l][k].conjugate();
		return C_mlk_conj;
	}

	ex fn_D (int m, int l, int k){
		ex D_mlk;
 /**
		 **  (q_l - q_k)^2
		 **     D_mlk = - 4 ---------------
		 **     AC_lk^2
  **/
		ex denom = mat_AC[l][k]*mat_AC[l][k];
		check0denom(denom, "D", m, l, k);

		ex temp = 0;
		temp = pow(mat_q[l][0] - mat_q[k][0], 2);
		for(int ii = 1; ii < 4; ii++){
			ex temp2 = mat_q[l][ii] - mat_q[k][ii];
			temp -= temp2*temp2;
		}
		D_mlk = -(temp / denom) * 4.0;
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
 /**    /-       d_lk
		 **     | = 0, if Im( - ------- ) < 0;
		 **     |AC_lk
		 **     |
		 **     |d_lk
		 **     f_lk    { = 1, if Im( - ------- ) = 0;
		 **     |AC_lk
		 **     |
		 **     |d_lk
		 **     | = 2, if Im( - ------- ) > 0;
		 **     |AC_lk
		 **     \-
  **/
		ex denom = mat_AC[l][k];
		check0denom(denom, "f", l, k);


		ex factor = - mat_d_im[l][k]/denom;
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
 /**    /-       d_lk
		 **     | = 0, if Im( - ------- ) > 0;
		 **     |AC_lk
		 **     |
		 **     |d_lk
		 **     fminus_lk       { = 1, if Im( - ------- ) = 0;
		 **     |AC_lk
		 **     |
		 **     |d_lk
		 **     | = 2, if Im( - ------- ) < 0;
		 **     |AC_lk
		 **     \-
  **/
		ex denom = mat_AC[l][k];
		check0denom(denom, "fminus", l, k);

		ex factor = - mat_d_im[l][k]/denom;
		return myfn_fminus(factor);
	}
	
	void xloopsGiNaC_calc_lev2(){
		int k, l, m;
#ifdef _AGK_debugmode
		printf("Level 2 Calculating. Please Wait a Long Moment ...\n");
#endif
		//	1.	Calc AC_lk & alpha_lk
		for(k = 0; k<4; k++) for(l = 0; l<4; l++){
			if (l!=k){ // calculate terms with dif. indices
#ifdef _AGK_debugmode
				printf("%d-%d\n", k+1, l+1);
#endif
				mat_AC[l][k] = fn_AC(l, k);
				mat_alpha[l][k] = fn_alpha(l, k);
			}
			else{ // set 0 value for same indices
				mat_AC[l][k] = 0;
				mat_alpha[l][k] = 0;
			}
		}
		
		//	2.	Calc f_lk, f_minus
		for(k = 0; k<4; k++) for(l = 0; l<4; l++){
			if (l!=k){ // calculate terms with dif. indices
#ifdef _AGK_debugmode
				printf("%d-%d\n", k+1, l+1);
#endif
				mat_f[l][k] = fn_f(l, k);
				mat_fminus[l][k] = fn_fminus(l,k);
			}
			else{ // set 0 value for same indices
				mat_f[l][k] = 0;
				mat_fminus[l][k] = 0;
			}
		}
		
		//	3.	Calc A_mlk, B_mlk, C_mlk (+ real, imag, conj), D_mlk
		for(k = 0; k<4; k++) for(l = 0; l<4; l++) for(m = 0; m<4; m++){
			if(m!=l && m!=k && l!=k){ // calculate terms with dif. indices
#ifdef _AGK_debugmode
				printf("%d-%d-%d\n", k+1, l+1, m+1);
#endif
				mat_A[m][l][k] = fn_A(m, l, k);
				mat_B[m][l][k] = fn_B(m, l, k);
				mat_C[m][l][k] = fn_C(m, l, k); mat_C_im[m][l][k] = fn_C_im(m, l, k); mat_C_re[m][l][k] = fn_C_re(m, l, k); mat_C_conj[m][l][k] = fn_C_conj(m, l, k);
				mat_D[m][l][k] = fn_D(m, l, k);
			}
			else{ // set 0 value for same indices
				mat_A[m][l][k] = 0;
				mat_B[m][l][k] = 0;
				mat_C[m][l][k] = 0; mat_C_im[m][l][k] = 0; mat_C_re[m][l][k] = 0; mat_C_conj[m][l][k] = 0;
				mat_D[m][l][k] = 0;
			}
		}
	}

}// Namespace xloops

