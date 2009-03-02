/*******************************************************************************
**
** xloop-ginacs Project
** 1Loop4Pt implementation
**
**
** HCMUNS,
** June 2008
**
**
** Author(s): Son, D. H. (sondo01@gmail.com)
** Khiem, Phan (phanhongkhiem@gmail.com)
** Quan, Phan (anhquan.phanle@gmail.com)
********************************************************************************
**
** level 3 = 2 dimensions D0 integral
**
********************************************************************************
**
** Historial Log:
** Date Version Author Description
** _________ _______ _________ ________________________________
** 200902211.0 Quan Phan Create this file
*******************************************************************************/

#include "ginac/ginac.h"
#include "extrm.h"
#include "my_fns.h"
#include "trmchk.h"
#include "lev3.h"

using namespace GiNaC;
namespace xloops{

	ex fn_beta (int m, int l, int k){
 /**A_mlk
		 ** term1 = ------- - alpha_lk
		 ** B_mlk
		 **
		 ** term1 + sqrt(term1^2 - D_mlk + i*rho)
		 ** beta_mlk = ----------------------------------------
		 ** D_mlk
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
 /**A_mlk
		 ** term1 = ------- - alpha_lk
		 ** B_mlk
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
 /** /- -C_mlk
		 ** | = 0, if Im( -------- ) < 0;
		 ** | B_mlk
		 ** |
		 ** | -C_mlk
		 ** g_mlk | = 1, if Im( -------- ) = 0;
		 ** | B_mlk
		 ** |
		 ** | -C_mlk
		 ** | = 2, if Im( -------- ) > 0;
		 ** | B_mlk
		 ** \-
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
 /** /- -C_mlk
		 ** | = 0, if Im( -------- ) > 0;
		 ** | B_mlk
		 ** |
		 ** | -C_mlk
		 ** gminus_mlk { = 1, if Im( -------- ) = 0;
		 ** | B_mlk
		 ** |
		 ** | -C_mlk
		 ** | = 2, if Im( -------- ) < 0;
		 ** | B_mlk
		 ** \-
  **/
		ex denom = mat_B[m][l][k];
		check0denom(denom, "gminus", m, l, k);

		ex factor = -mat_C_im[m][l][k]/denom;

		return myfn_gminus(factor);
	}

	ex fn_Q (int m, int l, int k){
		ex Q_mlk;
 /** C_mlk d_lk
		 ** Q_mlk = -2 ( ------- + ------- beta_mlk )
		 ** B_mlk AC_lk
  **/
		check0denom(mat_B[m][l][k], "Q", m, l, k);
		check0denom(mat_AC[l][k], "Q", m, l, k);
 
		Q_mlk = -( mat_C[m][l][k]/mat_B[m][l][k] + mat_beta[m][l][k]*mat_d[l][k]/ mat_AC[l][k] ) * 2.0;
		return Q_mlk;
	}
	ex fn_Q_im (int m, int l, int k){
		ex Q_mlk_im;
 
		Q_mlk_im = imag_part(mat_Q[m][l][k]);
		return Q_mlk_im;
	}
	ex fn_Q_re (int m, int l, int k){
		ex Q_mlk_re;
 
		Q_mlk_re = real_part(mat_Q[m][l][k]);
		return Q_mlk_re;
	}
	ex fn_Q_conj (int m, int l, int k){
		ex Q_mlk_conj;
 
		Q_mlk_conj = mat_Q[m][l][k].conjugate();
		return Q_mlk_conj;
	}

	ex fn_P (int m, int l, int k){
		ex P_mlk;
 /** A_mlk 
		 ** P_mlk = -2 [( ------- - alpha_lk) * (1 + beta_mlk.phi_mlk) - D_mlk . beta_mlk - phi_mlk ]
		 ** B_mlk
  **/
		check0denom(mat_B[m][l][k], "P", m, l, k);
 
		P_mlk = -(
			  (mat_A[m][l][k] / mat_B[m][l][k] - mat_alpha[l][k]) *(1 + mat_beta[m][l][k]*mat_phi[m][l][k])
				- mat_D[m][l][k]*mat_beta[m][l][k] 
				- mat_phi[m][l][k]
			 ) * 2.0;
		return P_mlk;
	}

	ex fn_E (int m, int l, int k){
		ex E_mlk;
 /** d_lk C_mlk
		 ** E_mlk = -2 ( ------- + ------- phi_mlk )
		 ** AC_lk B_mlk
  **/
		check0denom(mat_B[m][l][k], "E", m, l, k);
		check0denom(mat_AC[l][k], "E", m, l, k);
 
		E_mlk = -(
			  mat_d[l][k]/mat_AC[l][k] 
				+ mat_phi[m][l][k] * mat_C[m][l][k]/mat_B[m][l][k]
			 ) * 2.0;
		return E_mlk;
	}
	ex fn_E_im (int m, int l, int k){
		ex E_mlk_im;
 
		E_mlk_im = imag_part(mat_E[m][l][k]);
		return E_mlk_im;
	}
	ex fn_E_re (int m, int l, int k){
		ex E_mlk_re;
 
		E_mlk_re = real_part(mat_E[m][l][k]);
		return E_mlk_re;
	}

	ex fn_F (int n, int m, int l, int k){
		ex F_nmlk;
 /** C_nlk.B_mlk - C_mlk.B_nlk
		 ** F_nmlk = ---------------------------
		 ** A_nlk.B_mlk - A_mlk.B_nlk
  **/
		ex denom = mat_A[n][l][k] * mat_B[m][l][k] - mat_A[m][l][k] * mat_B[n][l][k];
		check0denom(denom, "F", m, l, k);
 
		F_nmlk = mat_C[n][l][k]*mat_B[m][l][k] - mat_C[m][l][k]*mat_B[n][l][k];
		F_nmlk /= denom;
		return F_nmlk;
	}
	ex fn_F_im (int n, int m, int l, int k){
		ex F_nmlk_im;
 
		F_nmlk_im = imag_part(mat_F[n][m][l][k]);
		return F_nmlk_im;
	}

	void xloopsGiNaC_calc_lev3(){
		// init and calculate all the terms. Use the same scheme: same index = set 0 value, dif index = equivalent function value
		int k, l, m, n;
		
#ifdef _AGK_debugmode
		printf("Level 3 Calculating. Please Wait a Long Moment ...\n");
#endif
		
		// calculate beta, phi, g, gminus
		for(k = 0; k<4; k++) for(l = 0; l<4; l++) for(m = 0; m<4; m++){
			if( m!= l && m!=k && l!=k){ // calc. dif indices
#ifdef _AGK_debugmode
				printf("%d-%d-%d\n", k+1, l+1, m+1);
#endif
				mat_beta[m][l][k] = fn_beta(m, l, k);
				mat_phi[m][l][k] = fn_phi(m, l, k);
				mat_g[m][l][k] = fn_g(m, l, k);
				mat_gminus[m][l][k] = fn_gminus(m, l, k);
			}
			else{ // set 0
				mat_beta[m][l][k] = 0;
				mat_phi[m][l][k] = 0;
				mat_g[m][l][k] = 0;
				mat_gminus[m][l][k] = 0;
			}
		}
		
		// calculate Q, P, E
		for(k=0; k<4; k++) for(l=0; l<4; l++) for(m = 0; m<4; m++){
			if(m!=l && m!=k && l!=k){
#ifdef _AGK_debugmode
				printf("%d-%d-%d\n", k+1, l+1, m+1);
#endif
				mat_Q[m][l][k]  =fn_Q(m, l, k); mat_Q_im[m][l][k]=fn_Q_im(m, l, k); mat_Q_re[m][l][k]=fn_Q_re(m, l, k); mat_Q_conj[m][l][k]=fn_Q_conj(m, l, k);
				mat_P[m][l][k]  =fn_P(m, l, k);
				mat_E[m][l][k]  =fn_E(m, l, k); mat_E_im[m][l][k]=fn_E_im(m, l, k); mat_E_re[m][l][k]=fn_E_re(m, l, k);
			}
			else{
				mat_Q[m][l][k] = 0; mat_Q_im[m][l][k] = 0; mat_Q_re[m][l][k] = 0; mat_Q_conj[m][l][k] = 0;
				mat_P[m][l][k] = 0;
				mat_E[m][l][k] = 0;
				mat_E_im[m][l][k] = 0; mat_E_re[m][l][k] = 0;
			}
		}
		
		// Calculate F
		for(k=0; k<4; k++) for(l=0; l<4; l++) for(m=0; m<4; m++) for(n=0; n<4; n++){
			if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k){
#ifdef _AGK_debugmode
				printf("%d-%d-%d-%d\n", k+1, l+1, m+1, n+1);
#endif
				mat_F[n][m][l][k] = fn_F(n, m, l, k); mat_F_im[n][m][l][k] = fn_F_im(n, m, l, k);
			}
			else{
				mat_F[n][m][l][k] = 0;
				mat_F_im[n][m][l][k] = 0;
			}
		}
	}
}// Namespace xloops
