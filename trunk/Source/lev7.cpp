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
**	Level 7 functions
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
#include "lev7.h"

using namespace GiNaC;
namespace xloops{
	//	II.7	Level 7 Variable Functions
	ex fn_OPlus (int n, int m, int l, int k){
		ex OPlus_nmlk;
	 /**                                               F_nmlk                     1 - beta_mlk.phi_mlk
		 ** OPlus_nmlk = (f_lk.g_mlk+fminus_lk.g_mlk).ln(----------) - f_lk.g_mlk.ln(----------------------) +
		 **                                               beta_mlk                           beta_mlk
		 **                           1 - beta_mlk.phi_mlk
		 **     - f_lk.gminus_mlk.ln(----------------------) - 2pi.i.f_lk.g_mlk.Sign(Im(F_nmlk)) +
		 **                               - beta_mlk
		 **                                           -P_mlk                                                 Q_mlk
		 **     - (f_lk.g_mlk + fminus_lk.g_mlk)[ ln(----------) - ln(P_mlk) + 2.i.pi.theta(-P_mlk).Sign(Im(-------))
		 **                                           beta_mlk                                               P_mlk
		 **                                                                P_mlk
		 **                                                + 2i.pi.theta(----------) + eta(z-z1beta_mlk,z-z2beta_mlk)] + 
		 **                                                               beta_mlk
		 **                                                                               Q_mlk
		 **     + f_lk.g_mlk[ln(-P_mlk.phi_mlk) - ln(P_mlk) + 2i.pi.theta(-P_mlk).Sign(Im-------) +
		 **                                                                               P_mlk
		 **                                               + 2i.pi.theta(P_mlk.phi_mlk) + eta(z-z1phi_mlk, z-z2phi_mlk)] + 
		 **                                                                                    Q_mlk
		 **     + f_lk.gminus_mlk[ln(-P_mlk.phi_mlk) - ln(-P_mlk) + 2i.pi.theta(P_mlk).Sign(Im-------) +
		 **                                                                                    P_mlk
		 **                                               + 2i.pi.theta(P_mlk.phi_mlk) + eta(z-z1phi_mlk, z-z2phi_mlk)]
	  **/
		check0denom(mat_beta[m][l][k], "OPlus", n, m,l ,k);
		check0denom(mat_P[m][l][k], "OPlus", n, m,l ,k);
	 
		OPlus_nmlk = (mat_f[l][k]*mat_g[m][l][k] + mat_fminus[l][k]*mat_g[m][l][k])*log((mat_F[n][m][l][k] + I*Rho*mat_beta[m][l][k])/mat_beta[m][l][k]);
		OPlus_nmlk += -(mat_f[l][k]*mat_g[m][l][k] + mat_fminus[l][k]*mat_g[m][l][k])*fn_eta_plus(-1.0/mat_beta[m][l][k], mat_z1beta[m][l][k], mat_z2beta[m][l][k], mat_P[m][l][k]);
		OPlus_nmlk += mat_f[l][k]*mat_g[m][l][k]*fn_eta_plus(-mat_phi[m][l][k], mat_z1phi[m][l][k], mat_z2phi[m][l][k], mat_P[m][l][k]);
		OPlus_nmlk += mat_f[l][k]*mat_gminus[m][l][k]*fn_eta_minus(-mat_phi[m][l][k], mat_z1phi[m][l][k], mat_z2phi[m][l][k], mat_P[m][l][k]);
	
		return OPlus_nmlk;
	}

	ex fn_OMinus (int n, int m, int l, int k){
		ex OMinus_nmlk;
	 /**                                         F_nmlk                           -F_nmlk
		 ** OMinus_nmlk = -fminus_lk.gminus_mlk.ln(----------) - f_lk.gminus_mlk.ln(----------) +
		 **                                         beta_mlk                         beta_mlk
	
		 **     - 2pi.i(fminus_lk.gminus_mlk+fminus_lk.g_mlk).Sign(Im(F_nmlk)) +
		 **                                                    1 - beta_mlk.phi_mlk
		 **     + (fminus_lk.gminus_mlk + fminus_lk.g_mlk).ln(----------------------) +
		 **                                                          beta_mlk
		 **                                   -P_mlk                                              Q_mlk
		 **     + fminus_lk.gminus_mlk.[ln(-----------) - ln(P_mlk) + 2i.pi.theta(-P_mlk).Sign(Im-------) +
		 **                                  beta_mlk                                             P_mlk
		 **                                                      P_mlk
		 **                                      + 2i.pi.theta(----------) + eta(z-z1beta_mlk, z-z2beta_mlk) ] +
		 **                                                     beta_mlk
		 **
		 **                             -P_mlk                                               Q_mlk
		 **     + f_lk.gminus_mlk.[ln(-----------) - ln(-P_mlk) + 2i.pi.theta(P_mlk).Sign(Im-------) +
		 **                            beta_mlk                                              P_mlk
		 **                                                      P_mlk
		 **                                      + 2i.pi.theta(----------) + eta(z-z1beta_mlk, z-z2beta_mlk) ] +
		 **                                                     beta_mlk
		 **
		 **                                                                                                              Q_mlk
		 **     + (fminus_lk.gminus_mlk + fminus_lk.g_mlk).[ln(-P_mlk.phi_mlk) - ln(P_mlk) + 2i.pi.theta(-P_mlk).Sign(Im-------) +
		 **                                                                                                              P_mlk
		 **
		 **                                      + 2i.pi.theta(P_mlk.phi_mlk) + eta(z-z1phi_mlk, z-z2phi_mlk) ]
		 **
	  **/
		check0denom(mat_beta[m][l][k], "OMinus", n, m,l ,k);
		check0denom(mat_P[m][l][k], "OMinus", n, m,l ,k);
		
		OMinus_nmlk = -mat_fminus[l][k]*mat_gminus[m][l][k]*log((mat_F[n][m][l][k] + I*Rho*mat_beta[m][l][k])/mat_beta[m][l][k])
				-mat_f[l][k]*mat_g[m][l][k]*log((mat_F[n][m][l][k]+I*Rho*mat_beta[m][l][k])/(-mat_beta[m][l][k]))
				+mat_fminus[l][k]*mat_gminus[m][l][k]*fn_eta_plus(-1.0/mat_beta[m][l][k], mat_z3beta[m][l][k], mat_z4beta[m][l][k], mat_P[m][l][k])
				+mat_f[l][k]*mat_gminus[m][l][k]*fn_eta_minus(-1.0/mat_beta[m][l][k], mat_z3beta[m][l][k], mat_z4beta[m][l][k], mat_P[m][l][k])
				-(mat_fminus[l][k]*mat_gminus[m][l][k] + mat_fminus[l][k]*mat_g[m][l][k])*fn_eta_plus(-mat_phi[m][l][k], mat_z3phi[m][l][k], mat_z4phi[m][l][k], mat_P[m][l][k]);
	
		return OMinus_nmlk;
	}

}// Namespace xloops
