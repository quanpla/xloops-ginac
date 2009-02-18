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
**	Level 3 functions
**
********************************************************************************
**
** Historial Log:
**	Date		Version	Author		Description
**	_________	_______	_________	________________________________
**	20090214	1.0	Quan Phan	Create this file
*******************************************************************************/

#include "ginac/ginac.h"
#include "extrm.h"
#include "my_fns.h"
#include "trmchk.h"
#include "lev3.h"

using namespace GiNaC;
namespace xloops{

	//	II.3	Level 3 Variable Functions
	ex fn_A (int m, int l, int k){
		ex A_mlk;
	 /**
		 **			     AC_mk
		 **	A_mlk = a_mk - a_lk -------
		 **			     AC_lk
	  **/
		ex denom = mat_AC[l][k];
		check0denom(denom, "A", m, l, k);
	
		A_mlk = mat_a[m][k] - mat_a[l][k]*mat_AC[m][k]/mat_AC[l][k];
		return A_mlk;
	}

	ex fn_B (int m, int l, int k){
		ex B_mlk;
	 /**
		 **			     AC_mk
		 **	B_mlk = b_mk - b_lk -------
		 **			     AC_lk
	  **/
		ex denom = mat_AC[l][k];
		check0denom(denom, "B", m, l, k);
	
		B_mlk = mat_b[m][k] - mat_b[l][k]*mat_AC[m][k]/mat_AC[l][k];
		return B_mlk;
	}

	ex fn_C (int m, int l, int k){
		ex C_mlk;
	 /**
		 **			     AC_mk
		 **	C_mlk = d_mk - d_lk -------
		 **			     AC_lk
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
		 **		     (q_l - q_k)^2
		 **	D_mlk = - 4 ---------------
		 **			AC_lk^2
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
	 /**		/-	         d_lk
		 **	       	| = 0, if Im( - ------- ) < 0;
		 **		|		 AC_lk
		 **		|
		 **	       	|	         d_lk
		 **	f_lk	{ = 1, if Im( - ------- ) = 0;
		 **	  	|		 AC_lk
		 **		|
		 **		|	         d_lk
		 **		| = 2, if Im( - ------- ) > 0;
		 **		|		 AC_lk
		 **		\-
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
	 /**			/-	         d_lk
		 **	       		| = 0, if Im( - ------- ) > 0;
		 **			|		 AC_lk
		 **			|
		 **	 	      	|	         d_lk
		 **	fminus_lk	{ = 1, if Im( - ------- ) = 0;
		 **	  		|		 AC_lk
		 **			|
		 **			|	         d_lk
		 **			| = 2, if Im( - ------- ) < 0;
		 **			|		 AC_lk
		 **			\-
	  **/
		ex denom = mat_AC[l][k];
		check0denom(denom, "fminus", l, k);
		
		ex factor = - mat_d_im[l][k]/denom;
		return myfn_fminus(factor);
	}

}// Namespace xloops
