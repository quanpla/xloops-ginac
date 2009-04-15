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
#include "lev1.h"
#include "lev2.h"
#include "lev3.h"
#include "lev4.h"

using namespace GiNaC;
namespace xloops{
	ex fn_z1phi (int m, int l, int k){
		ex z1phi_mlk;
	// init
	   	ex E_mlk = fn_E(m, l, k), phi_mlk = fn_phi(m, l, k), Q_mlk = fn_Q(m, l, k), msquare_k = mat_msquare[k];

		ex denom = -2.0 * P_mlk * phi_mlk;

		check0denom(denom, "z1phi", m, l, k);

	// calculation
		z1phi_mlk = - E_mlk + Q_mlk*phi_mlk;
		z1phi_mlk += sqrt( pow(z1phi_mlk, 2) 
		- 4.0*P_mlk*phi_mlk*(msquare_k - I*Rho /*Feynman's pres.*/)
		);
		z1phi_mlk /= denom;

		return z1phi_mlk;
	}

	ex fn_z2phi (int m, int l, int k){
		ex z2phi_mlk;
	// init
	   	ex E_mlk = fn_E(m, l, k), phi_mlk = fn_phi(m, l, k), Q_mlk = fn_Q(m, l, k), msquare_k = mat_msquare[k];

		ex denom =  -2.0 * P_mlk * phi_mlk;

		check0denom(denom, "z2phi", m, l, k);

	// calculation
		z2phi_mlk =  - E_mlk + Q_mlk*phi_mlk;
		z2phi_mlk -= sqrt(
				  pow(z2phi_mlk, 2)
				- 4.0*P_mlk*phi_mlk*(msquare_k - I*Rho /*Feynman's pres.*/)
				 );
		z2phi_mlk /= denom;

		return z2phi_mlk;
	}

	ex fn_z3phi (int m, int l, int k){
		ex z3phi_mlk;
	// init
	   	ex E_mlk = fn_E(m, l, k), phi_mlk = fn_phi(m, l, k), Q_mlk = fn_Q(m, l, k), msquare_k = mat_msquare[k];

		ex denom = -2.0 * P_mlk * phi_mlk;

		check0denom(denom, "z3phi", m, l, k);

	// calculation
		z3phi_mlk = E_mlk - Q_mlk*phi_mlk;
		z3phi_mlk += sqrt(
				  pow(z3phi_mlk, 2)
				- 4.0*P_mlk*phi_mlk*(msquare_k - I*Rho /*Feynman's pres.*/)
				 );
		z3phi_mlk /= denom;

		return z3phi_mlk;
	}

	ex fn_z4phi (int m, int l, int k){
		ex z4phi_mlk;
		// init
	   	ex E_mlk = fn_E(m, l, k), phi_mlk = fn_phi(m, l, k), Q_mlk = fn_Q(m, l, k), msquare_k = mat_msquare[k];

		ex denom =  -2.0 * P_mlk * phi_mlk;

		check0denom(denom, "z4phi", m, l, k);

// calculation
		z4phi_mlk = E_mlk - Q_mlk*phi_mlk;
		z4phi_mlk -= sqrt(
				  pow(z4phi_mlk, 2)
				- 4.0*P_mlk*phi_mlk*(msquare_k - I*Rho /*Feynman's pres.*/)
				 );
		z4phi_mlk /= denom;

		return z4phi_mlk;
	}

	ex fn_z1beta (int m, int l, int k){
		ex z1beta_mlk;
		// init
	   	ex E_mlk = fn_E(m, l, k), phi_mlk = fn_phi(m, l, k), Q_mlk = fn_Q(m, l, k), msquare_k = mat_msquare[k];

		check0denom(beta_mlk, "z1beta", m, l, k);
		check0denom(P_mlk, "z1beta", m, l, k);

// calculation
		z1beta_mlk = - E_mlk + Q_mlk/beta_mlk;
		z1beta_mlk += sqrt(
				   pow(z1beta_mlk, 2)
				- 4.0*P_mlk*(msquare_k - I*Rho/*Feynman's pres.*/)/beta_mlk
				  );
		z1beta_mlk /= - 2.0 * P_mlk / beta_mlk;
		return z1beta_mlk;
	}
	ex fn_z2beta (int m, int l, int k){
		ex z2beta_mlk;
		// init
	   	ex E_mlk = fn_E(m, l, k), phi_mlk = fn_phi(m, l, k), Q_mlk = fn_Q(m, l, k), msquare_k = mat_msquare[k];

		check0denom(beta_mlk, "z2beta", m, l, k);
		check0denom(P_mlk, "z2beta", m, l, k);

		// calculation
		z2beta_mlk = - E_mlk + Q_mlk/beta_mlk;
		z2beta_mlk -= sqrt(
				   pow(z2beta_mlk, 2)
				- 4.0*P_mlk*(msquare_k - I*Rho/*Feynman's pres.*/)/beta_mlk
				  );
		z2beta_mlk /= - 2.0 * P_mlk / beta_mlk;
		return z2beta_mlk;
	}


	ex fn_z3beta (int m, int l, int k){
		ex z3beta_mlk;
	   	ex E_mlk = fn_E(m, l, k), phi_mlk = fn_phi(m, l, k), Q_mlk = fn_Q(m, l, k), msquare_k = mat_msquare[k];

		check0denom(beta_mlk, "z3beta", m, l, k);
		check0denom(P_mlk, "z3beta", m, l, k);

		z3beta_mlk = E_mlk - Q_mlk/beta_mlk;
		z3beta_mlk += sqrt(
				   pow(z3beta_mlk, 2)
				- 4.0*P_mlk*(msquare_k - I*Rho/*Feynman's pres.*/)/beta_mlk
				  );
		z3beta_mlk /= - 2.0 * P_mlk / beta_mlk;
		return z3beta_mlk;
	}

	ex fn_z4beta (int m, int l, int k){
		ex z4beta_mlk;
	   	ex E_mlk = fn_E(m, l, k), phi_mlk = fn_phi(m, l, k), Q_mlk = fn_Q(m, l, k), msquare_k = mat_msquare[k];

		check0denom(beta_mlk, "z4beta", m, l, k);
		check0denom(P_mlk, "z4beta", m, l, k);

		z4beta_mlk = E_mlk - Q_mlk/beta_mlk;
		z4beta_mlk -= sqrt(
				   pow(z4beta_mlk, 2)
				- 4.0*P_mlk*(msquare_k - I*Rho/*Feynman's pres.*/)/beta_mlk
				  );
		z4beta_mlk /= - 2.0 * P_mlk / beta_mlk;
		return z4beta_mlk;
	}


	ex fn_T1 (int n, int m, int l, int k){
		ex T1_nmlk;

	   ex P_mlk = fn_P(m, l, k), Q_mlk = fn_Q(m, l, k), F_nmlk = fn_F(n, m, l, k), beta_mlk = fn_beta(m, l, k), E_mlk = fn_E(m, l, k), msquare_k = mat_msquare[k];
		check0denom(P_mlk, "T1", n, m, l, k);

		T1_nmlk = Q_mlk + P_mlk*F_nmlk - beta_mlk*E_mlk;
		T1_nmlk += sqrt(
				pow(T1_nmlk,2)
				- 4.0 * P_mlk*(
						Q_mlk*F_nmlk
						+ beta_mlk*msquare_k
						- I*beta_mlk*Rho /*Feynman's pres.*/
						       )
			       );
		T1_nmlk /= -2.0 * P_mlk;

		return T1_nmlk;
	}

	ex fn_T2 (int n, int m, int l, int k){
		ex T2_nmlk;
		check0denom(P_mlk, "T2", n, m, l, k);

		T2_nmlk = Q_mlk + P_mlk*F_nmlk - beta_mlk*E_mlk;
		T2_nmlk -= sqrt(
				T2_nmlk*T2_nmlk
				- 4.0 * P_mlk * (
				Q_mlk*F_nmlk
				+ beta_mlk*msquare_k
				- I*beta_mlk*Rho /*Feynman's pres.*/
							 )
			       );
		T2_nmlk /= -2.0 * P_mlk;
		return T2_nmlk;
	}
	ex fn_T3 (int n, int m, int l, int k){
		ex T3_nmlk;
		check0denom(P_mlk, "T3", n, m, l, k);

		T3_nmlk = -(Q_mlk + P_mlk*F_nmlk - beta_mlk*E_mlk);
		T3_nmlk += sqrt(
				pow(T3_nmlk,2)
				- 4.0 * P_mlk*(
						Q_mlk*F_nmlk
						+ beta_mlk*msquare_k
						- I*beta_mlk*Rho /*Feynman's pres.*/
						       )
			       );
		T3_nmlk /= -2.0 * P_mlk;

		return T3_nmlk;
	}

	ex fn_T4 (int n, int m, int l, int k){
		ex T4_nmlk;
		check0denom(P_mlk, "T4", n, m, l, k);

		T4_nmlk = - ( Q_mlk + P_mlk*F_nmlk - beta_mlk*E_mlk);
		T4_nmlk -= sqrt(
				pow(T4_nmlk,2)
				- 4.0 * P_mlk * (
						Q_mlk*F_nmlk
						+ beta_mlk*msquare_k
						- I*beta_mlk*Rho /*Feynman's pres.*/
							 )
			       );
		T4_nmlk /= -2.0 * P_mlk;
		return T4_nmlk;
	} 
}// Namespace xloops
