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
**	20090413	1.0	Quan Phan	can get value directly
**	20090214	1.0	Quan Phan	Create this file
*******************************************************************************/

#include "ginac/ginac.h"
#include "extrm.h"
#include "my_fns.h"
#include "trmchk.h"
#include "lev1.h"
#include "lev2.h"
#include "lev3.h"

using namespace GiNaC;
namespace xloops{

	ex fn_beta (int m, int l, int k){
		ex beta_mlk;
		
		ex A_mlk = fn_A(m, l, k), B_mlk = fn_B(m, l, k), D_mlk = fn_D(m, l, k), alpha_lk = fn_alpha(l, k);

		check0denom(D_mlk, "beta", m, l, k);
		check0denom(B_mlk, "beta", m, l, k);
	
		ex term1 = A_mlk / B_mlk - alpha_lk;
	
		beta_mlk = term1 + sqrt(term1*term1 - D_mlk);
		beta_mlk = beta_mlk / D_mlk;
	
		return beta_mlk;
	}

	ex fn_phi (int m, int l, int k){
		ex phi_mlk;
		ex A_mlk = fn_A(m, l, k), B_mlk = fn_B(m, l, k), D_mlk = fn_D(m, l, k), alpha_lk = fn_alpha(l, k);
		
		check0denom(B_mlk, "phi", m, l, k);
	
		ex term1 = A_mlk/B_mlk - alpha_lk;
		phi_mlk = term1 + sqrt(term1*term1 - D_mlk);
	
		return phi_mlk;
	}

///TODO: Remove if not use anymore
	ex myfn_g_eval(const ex &factor){
		if (is_a<numeric>(factor)){
			if (my_is_zero(factor) == 1.0) // if zero, return 1
				return 1.0;
			if (factor.info(info_flags::negative)) // if negative, return 0
				return 0.0;
			return 2.0; // if positive return 2
		}
		else{
			// if not a number, do not calculate
			return myfn_g(factor).hold();
		}
	}

	REGISTER_FUNCTION(myfn_g, eval_func(myfn_g_eval));
	ex fn_g (int m, int l, int k){
		
		ex B_mlk = fn_B(m, l, k), C_im_mlk = fn_C_im(m, l, k);

		check0denom(B_mlk, "g", m, l, k);

		ex factor = -C_im_mlk/B_mlk;
		return myfn_g(factor);
	}

	ex myfn_gminus_eval(const ex &factor){
		if (is_a<numeric>(factor)){
			if (my_is_zero(factor) == 1.0) // if zero, return 1
				return 1.0;
			if (factor.info(info_flags::negative)) // if negative, return 2
				return 2.0;
			return 0.0; // if positive, return 0
		}
		else{
			// if not a number, don't calculate
			return myfn_gminus(factor).hold();
		}
	}

	REGISTER_FUNCTION(myfn_gminus, eval_func(myfn_gminus_eval));
	ex fn_gminus (int m, int l, int k){
	  	ex B_mlk = fn_B(m, l, k), C_im_mlk = imag_part( fn_C(m, l, k) );
	  	
		check0denom(B_mlk, "gminus", m, l, k);

		ex factor = -C_im_mlk/B_mlk;
		return myfn_gminus(factor); // read myfn_gminus for the implementation of this function
	}
 
	ex fn_Q (int m, int l, int k){
		ex Q_mlk;
 /** C_mlk d_lk
		 ** Q_mlk = -2 ( ------- + ------- beta_mlk )
		 ** B_mlk AC_lk
  **/
  		ex B_mlk = fn_B(m, l, k), AC_lk = fn_AC(l, k), C_mlk = fn_C(m, l, k), beta_mlk = fn_beta(m, l, k), d_lk = fn_d(l, k);
  		
  		check0denom(B_mlk, "Q", m, l, k);
		check0denom(AC_lk, "Q", m, l, k);
 
		Q_mlk = -2.0*( C_mlk/B_mlk + beta_mlk*d_lk/AC_lk );
		return Q_mlk;
	}
	ex fn_Q_im (int m, int l, int k){
		//init
		ex Q_mlk_im;
 		ex Q_mlk = fn_Q(m, l, k);

		//calc.
		Q_mlk_im = imag_part(Q_mlk);
		return Q_mlk_im;
	}
	ex fn_Q_re (int m, int l, int k){
		//init
		ex Q_mlk_re;
 		ex Q_mlk = fn_Q(m, l, k);

		//calc.
		Q_mlk_re = real_part(Q_mlk);
		return Q_mlk_re;
	}
	ex fn_Q_conj (int m, int l, int k){
		ex Q_mlk_conj;
  		ex Q_mlk = fn_Q(m, l, k);

		//calc.
		Q_mlk_conj = Q_mlk.conjugate();
		return Q_mlk_conj;
	}

	ex fn_E (int m, int l, int k){
		ex E_mlk;
	  	ex B_mlk = fn_B(m, l, k), C_mlk = fn_C(m, l, k), AC_lk = fn_AC(l, k), d_lk = fn_d(l, k), phi_mlk = fn_phi(m, l, k);
		check0denom(B_mlk, "E", m, l, k);
		check0denom(AC_lk, "E", m, l, k);
 
		E_mlk = -2.0 * ( C_mlk*phi_mlk/B_mlk + d_lk/AC_lk );
		
		return E_mlk;
	}
	ex fn_E_im (int m, int l, int k){
		//init
		ex E_mlk_im;
 		ex E_mlk = fn_E(m, l, k);

		// calc 
		E_mlk_im = imag_part(E_mlk);
		return E_mlk_im;
	}
	ex fn_E_re (int m, int l, int k){
		// init
		ex E_mlk_re;
 		ex E_mlk = fn_E(m, l, k);

		// calc
		E_mlk_re = real_part(E_mlk);
		return E_mlk_re;
	}

	ex fn_P (int m, int l, int k){
		ex P_mlk;
  		ex A_mlk = fn_A(m, l, k), B_mlk = fn_B(m, l, k), D_mlk = fn_D(m, l, k), phi_mlk = fn_phi(m, l, k), alpha_lk = fn_alpha(l, k), beta_mlk = fn_beta(m, l, k);
  		
		check0denom(B_mlk, "P", m, l, k);
 
		P_mlk = -2.0 * (
			(A_mlk/B_mlk - alpha_lk) * (1.0 + beta_mlk*phi_mlk)
			- D_mlk*beta_mlk
			- phi_mlk
			);
		
		return P_mlk;
	}

	ex fn_A0 (int m, int l, int k){
		ex A0_mlk;
		ex P_mlk = fn_P(m, l, k), E_mlk = fn_E(m, l, k);
		
		// calc.
		A0_mlk = imag_part(P_mlk*E_mlk);
		return A0_mlk;
	}

	ex fn_B0 (int m, int l, int k){
		ex B0_mlk;
		ex Q_mlk_conj = fn_Q_conj(m, l, k), E_mlk = fn_E(m, l, k)
			msquare_k = mat_msquare[k], P_mlk = fn_P(m, l, k);
		
		// calc.
		B0_mlk = imag_part( Q_mlk_conj*E_mlk - msquare_k*P_mlk + I*Rho*P_mlk );
		return B0_mlk;
	}

	ex fn_C0 (int m, int l, int k){
		ex B0_mlk;
		ex Q_mlk_conj = fn_Q_conj(m, l, k), E_mlk = fn_E(m, l, k)
			msquare_k = mat_msquare[k], P_mlk = fn_P(m, l, k);
		
		// calc.
		B0_mlk = imag_part( -Q_mlk_conj*msquare_k + Q_mlk_conj*I*Rho);
		return B0_mlk;
	}

	ex fn_F (int n, int m, int l, int k){
		ex F_nmlk;
 /** C_nlk.B_mlk - C_mlk.B_nlk
		 ** F_nmlk = ---------------------------
		 ** A_nlk.B_mlk - A_mlk.B_nlk
  **/
  		ex A_nlk = fn_A(n, l, k), A_mlk = fn_A(m, l, k), B_nlk = fn_B(n, l, k), B_mlk = fn_B(m, l, k), C_nlk = fn_C(n, l, k), C_mlk = fn_C(m, l, k)
  			, beta_mlk = fn_beta(m, l, k);

		ex denom = A_nlk * B_mlk - A_mlk * B_nlk + I*Rho*beta_mlk;
		check0denom(denom, "F", m, l, k);
 
		F_nmlk = C_nlk*B_mlk - C_mlk*B_nlk;
		F_nmlk /= denom;
		return F_nmlk;
	}
	ex fn_F_im (int n, int m, int l, int k){
		// init
		ex F_nmlk_im;
 		ex F_nmlk = fn_F(n, m, l, k);

		// calc.
		F_nmlk_im = imag_part(F_nmlk);
		return F_nmlk_im;
	}

	
}// Namespace xloops
