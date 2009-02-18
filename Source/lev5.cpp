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
**	Level 5 functions
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
#include "lev5.h"

using namespace GiNaC;
namespace xloops{
	//	II.5	Level 5 Variable Functions
	ex fn_Q (int m, int l, int k){
		ex Q_mlk;
	 /**               C_mlk     d_lk
		 ** Q_mlk = -2 ( ------- + ------- beta_mlk )
		 **               B_mlk     AC_lk
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
	 	/**               A_mlk 
		 ** P_mlk = -2 [( ------- - alpha_lk) * (1 + beta_mlk.phi_mlk) - D_mlk . beta_mlk - phi_mlk ]
		 **                B_mlk
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
	 	/**               d_lk      C_mlk
		 ** E_mlk = -2 ( ------- + ------- phi_mlk )
		 **               AC_lk     B_mlk
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
	 	/**           C_nlk.B_mlk - C_mlk.B_nlk
		 ** F_nmlk = ---------------------------
		 **           A_nlk.B_mlk - A_mlk.B_nlk
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

}// Namespace xloops
