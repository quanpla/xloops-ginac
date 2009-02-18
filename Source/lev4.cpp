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
**	Level 4 functions
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
#include "lev4.h"

using namespace GiNaC;
namespace xloops{

	//	II.4	Level 4 Variable Functions
	ex fn_beta (int m, int l, int k){
	 /**	    A_mlk
		 ** term1 = ------- - alpha_lk
		 **	    B_mlk
		 **
		 **		term1 + sqrt(term1^2 - D_mlk + i*rho)
		 ** beta_mlk = ----------------------------------------
		 **				D_mlk
	  **/
		ex beta_mlk;
		ex denom = mat_D[m][l][k];
		check0denom(mat_D[m][l][k], "beta", m, l, k);
		check0denom(mat_B[m][l][k], "beta", m, l, k);
	
		ex term1 = mat_A[m][l][k] / mat_B[m][l][k] - mat_alpha[l][k];
	
		beta_mlk = term1 + sqrt(term1*term1 - mat_D[m][l][k] + I * Rho1 /*This is NOT Feynman's prescription.*/);
		beta_mlk = beta_mlk / mat_D[m][l][k];
	
		return beta_mlk;
	}

	ex fn_phi (int m, int l, int k){
	 /**	    A_mlk
		 ** term1 = ------- - alpha_lk
		 **	    B_mlk
		 **
		 ** beta_mlk = term1 + sqrt(term1^2 - D_mlk + i*eta)
	  **/
		ex phi_mlk;
		ex denom = mat_B[m][l][k];
		check0denom(denom, "phi", m, l, k);
	
		ex term1 = mat_A[m][l][k] / mat_B[m][l][k] - mat_alpha[l][k];
	
		phi_mlk = term1 + sqrt(term1*term1 - mat_D[m][l][k] + I * Rho1 /*This is NOT Feynman's prescription*/);
	
		return phi_mlk;
	}

	ex myfn_g_eval(const ex &factor){
		if (is_a<numeric>(factor)){
			if (my_is_zero(factor) == 1.0)
				return 1.0;
			if (factor.info(info_flags::negative))
				return 0.0;
			return 2.0;
		}
		else{
			return myfn_g(factor).hold();
		}
	}

	REGISTER_FUNCTION(myfn_g, eval_func(myfn_g_eval));
	ex fn_g (int m, int l, int k){
	 /**			/-	       -C_mlk
		 **	       		| = 0, if Im( -------- ) < 0;
		 **			|		B_mlk
		 **			|
		 **	 	      	|	       -C_mlk
		 **	g_mlk		| = 1, if Im( -------- ) = 0;
		 **	  		|		B_mlk
		 **			|
		 **			|	       -C_mlk
		 **			| = 2, if Im( -------- ) > 0;
		 **			|		B_mlk
		 **			\-
	  **/
		ex denom = mat_B[m][l][k];
		check0denom(denom, "g", m, l, k);

		ex factor = -mat_C_im[m][l][k]/denom;
		return myfn_g(factor);
	}

	ex myfn_gminus_eval(const ex &factor){
		if (is_a<numeric>(factor)){
			if (my_is_zero(factor) == 1.0)
				return 1.0;
			if (factor.info(info_flags::negative))
				return 2.0;
			return 0.0;
		}
		else{
			return myfn_gminus(factor).hold();
		}
	}

	REGISTER_FUNCTION(myfn_gminus, eval_func(myfn_gminus_eval));
	ex fn_gminus (int m, int l, int k){
	 /**			/-	       -C_mlk
		 **	       		| = 0, if Im( -------- ) > 0;
		 **			|		B_mlk
		 **			|
		 **	 	      	|	       -C_mlk
		 **	gminus_mlk	{ = 1, if Im( -------- ) = 0;
		 **	  		|		B_mlk
		 **			|
		 **			|	       -C_mlk
		 **			| = 2, if Im( -------- ) < 0;
		 **			|		B_mlk
		 **			\-
	  **/
		ex denom = mat_B[m][l][k];
		check0denom(denom, "gminus", m, l, k);

		ex factor = -mat_C_im[m][l][k]/denom;
	
		return myfn_gminus(factor);
	}

}// Namespace xloops
