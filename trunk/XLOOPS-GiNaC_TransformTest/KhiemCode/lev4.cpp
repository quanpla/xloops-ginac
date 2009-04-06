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
**Khiem, Phan     (phanhongkhiem@gmail.com)
**Quan, Phan      (anhquan.phanle@gmail.com)
********************************************************************************
**
**      Level 4 = 1 dimension D0 integral
**
********************************************************************************
**
** Historial Log:
**      Date     Version Author   Description
**      ________________ _________________________________________
**      20090221 1.0     Quan PhanCreate this file
*******************************************************************************/

#include "ginac/ginac.h"
#include "extrm.h"
#include "my_fns.h"
#include "trmchk.h"
#include "lev4.h"

using namespace GiNaC;
namespace xloops{

	ex fn_z1phi (int m, int l, int k){
		ex z1phi_mlk;
  /**
   **/
		ex denom = -2.0 * mat_P[m][l][k] * mat_phi[m][l][k];
 
		check0denom(denom, "z1phi", m, l, k);

		z1phi_mlk = - mat_E[m][l][k] + mat_Q[m][l][k]*mat_phi[m][l][k];
		z1phi_mlk += sqrt(
				  pow(z1phi_mlk, 2)
				- 4.0*mat_P[m][l][k]*mat_phi[m][l][k]*(mat_msquare[k] - I*Rho /*This is Feynman's prescription*/)
				 );
		z1phi_mlk /= denom;

		return z1phi_mlk;
	}

	ex fn_z2phi (int m, int l, int k){
		ex z2phi_mlk;
  /**
   **/
		ex denom =  -2.0 * mat_P[m][l][k] * mat_phi[m][l][k];

		check0denom(denom, "z2phi", m, l, k);

		z2phi_mlk =  - mat_E[m][l][k] + mat_Q[m][l][k]*mat_phi[m][l][k];
		z2phi_mlk -= sqrt(
				  pow(z2phi_mlk, 2)
				- 4.0*mat_P[m][l][k]*mat_phi[m][l][k]*(mat_msquare[k] - I*Rho /*This is Feynman's prescription*/)
				 );
		z2phi_mlk /= denom;

		return z2phi_mlk;
	}

	ex fn_z3phi (int m, int l, int k){
		ex z3phi_mlk;
  /**
   **/
		ex denom = -2.0 * mat_P[m][l][k] * mat_phi[m][l][k];

		check0denom(denom, "z3phi", m, l, k);

		z3phi_mlk = mat_E[m][l][k] - mat_Q[m][l][k]*mat_phi[m][l][k];
		z3phi_mlk += sqrt(
				  pow(z3phi_mlk, 2)
				- 4.0*mat_P[m][l][k]*mat_phi[m][l][k]*(mat_msquare[k] - I*Rho /*This is Feynman's prescription*/)
				 );
		z3phi_mlk /= denom;

		return z3phi_mlk;
	}

	ex fn_z4phi (int m, int l, int k){
		ex z4phi_mlk;
  /**
   **/
		ex denom =  -2.0 * mat_P[m][l][k] * mat_phi[m][l][k];

		check0denom(denom, "z4phi", m, l, k);

		z4phi_mlk = mat_E[m][l][k] - mat_Q[m][l][k]*mat_phi[m][l][k];
		z4phi_mlk -= sqrt(
				  pow(z4phi_mlk, 2)
				- 4.0*mat_P[m][l][k]*mat_phi[m][l][k]*(mat_msquare[k] - I*Rho /*This is Feynman's prescription*/)
				 );
		z4phi_mlk /= denom;

		return z4phi_mlk;
	}

	ex fn_z1beta (int m, int l, int k){
		ex z1beta_mlk;
 /**
  **/
		check0denom(mat_beta[m][l][k], "z1beta", m, l, k);
		check0denom(mat_P[m][l][k], "z1beta", m, l, k);

		z1beta_mlk = - mat_E[m][l][k] + mat_Q[m][l][k]/mat_beta[m][l][k];
		z1beta_mlk += sqrt(
				   pow(z1beta_mlk, 2)
				- 4.0*mat_P[m][l][k]*(mat_msquare[k] - I*Rho/*This is Feynman's prescription*/)/mat_beta[m][l][k]
				  );
		z1beta_mlk /= - 2.0 * mat_P[m][l][k] / mat_beta[m][l][k];
		return z1beta_mlk;
	}
	ex fn_z2beta (int m, int l, int k){
		ex z2beta_mlk;
 /**
  **/
		check0denom(mat_beta[m][l][k], "z2beta", m, l, k);
		check0denom(mat_P[m][l][k], "z2beta", m, l, k);

		z2beta_mlk = - mat_E[m][l][k] + mat_Q[m][l][k]/mat_beta[m][l][k];
		z2beta_mlk -= sqrt(
				   pow(z2beta_mlk, 2)
				- 4.0*mat_P[m][l][k]*(mat_msquare[k] - I*Rho/*This is Feynman's prescription*/)/mat_beta[m][l][k]
				  );
		z2beta_mlk /= - 2.0 * mat_P[m][l][k] / mat_beta[m][l][k];
		return z2beta_mlk;
	}


	ex fn_z3beta (int m, int l, int k){
		ex z3beta_mlk;
 /**
  **/
		check0denom(mat_beta[m][l][k], "z3beta", m, l, k);
		check0denom(mat_P[m][l][k], "z3beta", m, l, k);

		z3beta_mlk = mat_E[m][l][k] - mat_Q[m][l][k]/mat_beta[m][l][k];
		z3beta_mlk += sqrt(
				   pow(z3beta_mlk, 2)
				- 4.0*mat_P[m][l][k]*(mat_msquare[k] - I*Rho/*This is Feynman's prescription*/)/mat_beta[m][l][k]
				  );
		z3beta_mlk /= - 2.0 * mat_P[m][l][k] / mat_beta[m][l][k];
		return z3beta_mlk;
	}

	ex fn_z4beta (int m, int l, int k){
		ex z4beta_mlk;
 /**
  **/
		check0denom(mat_beta[m][l][k], "z4beta", m, l, k);
		check0denom(mat_P[m][l][k], "z4beta", m, l, k);

		z4beta_mlk = mat_E[m][l][k] - mat_Q[m][l][k]/mat_beta[m][l][k];
		z4beta_mlk -= sqrt(
				   pow(z4beta_mlk, 2)
				- 4.0*mat_P[m][l][k]*(mat_msquare[k] - I*Rho/*This is Feynman's prescription*/)/mat_beta[m][l][k]
				  );
		z4beta_mlk /= - 2.0 * mat_P[m][l][k] / mat_beta[m][l][k];
		return z4beta_mlk;
	}


	ex fn_T1 (int n, int m, int l, int k){
		ex T1_nmlk;
  /**     Q_mlk + P_mlk.F_nmlk - beta_mlk.E_mlk)
		 ** T1_nmlk = ---------------------------------------- +
		 **  - 2.P_mlk
		 ** sqrt((Q_mlk+P_mlk.F_nmlk-beta_mlk.E_mlk)^2 - 4.P_mlk(Q_mlk.F_nmlk+beta_mlk.m_k^2-i.beta_mlk.rho))
		 **      + -------------------------------------------------------------------------------------------------
		 **- 2.P_mlk
   **/
		check0denom(mat_P[m][l][k], "T1", n, m, l, k);

		T1_nmlk = mat_Q[m][l][k] + mat_P[m][l][k]*mat_F[n][m][l][k] - mat_beta[m][l][k]*mat_E[m][l][k];
		T1_nmlk += sqrt(
				pow(T1_nmlk,2)
				- 4.0 * mat_P[m][l][k]*(
						mat_Q[m][l][k]*mat_F[n][m][l][k]
						+ mat_beta[m][l][k]*mat_msquare[k]
						- I*mat_beta[m][l][k]*Rho /*This is Feynman's prescription*/
						       )
			       );
		T1_nmlk /= -2.0 * mat_P[m][l][k];

		return T1_nmlk;
	}

	ex fn_T2 (int n, int m, int l, int k){
		ex T2_nmlk;
  /**     Q_mlk + P_mlk.F_nmlk - beta_mlk.E_mlk)
		 ** T2_nmlk = ---------------------------------------- +
		 **  - 2.P_mlk
		 ** sqrt((Q_mlk+P_mlk.F_nmlk-beta_mlk.E_mlk)^2 - 4.P_mlk(Q_mlk.F_nmlk+beta_mlk.m_k^2-i.beta_mlk.rho))
		 **      - -------------------------------------------------------------------------------------------------
		 **- 2.P_mlk
   **/
		check0denom(mat_P[m][l][k], "T2", n, m, l, k);

		T2_nmlk = mat_Q[m][l][k] + mat_P[m][l][k]*mat_F[n][m][l][k] - mat_beta[m][l][k]*mat_E[m][l][k];
		T2_nmlk -= sqrt(
				T2_nmlk*T2_nmlk
				- 4.0 * mat_P[m][l][k] * (
				mat_Q[m][l][k]*mat_F[n][m][l][k]
				+ mat_beta[m][l][k]*mat_msquare[k]
				- I*mat_beta[m][l][k]*Rho /*This is Feynman's prescription*/
							 )
			       );
		T2_nmlk /= -2.0 * mat_P[m][l][k];
		return T2_nmlk;
	}
	ex fn_T3 (int n, int m, int l, int k){
		ex T3_nmlk;
  /**     Q_mlk + P_mlk.F_nmlk - beta_mlk.E_mlk)
		 ** T1_nmlk = ---------------------------------------- +
		 **  - 2.P_mlk
		 ** sqrt((Q_mlk+P_mlk.F_nmlk-beta_mlk.E_mlk)^2 - 4.P_mlk(Q_mlk.F_nmlk+beta_mlk.m_k^2-i.beta_mlk.rho))
		 **      + -------------------------------------------------------------------------------------------------
		 **- 2.P_mlk
   **/
		check0denom(mat_P[m][l][k], "T3", n, m, l, k);

		T3_nmlk = -(mat_Q[m][l][k] + mat_P[m][l][k]*mat_F[n][m][l][k] - mat_beta[m][l][k]*mat_E[m][l][k]);
		T3_nmlk += sqrt(
				pow(T3_nmlk,2)
				- 4.0 * mat_P[m][l][k]*(
						mat_Q[m][l][k]*mat_F[n][m][l][k]
						+ mat_beta[m][l][k]*mat_msquare[k]
						- I*mat_beta[m][l][k]*Rho /*This is Feynman's prescription*/
						       )
			       );
		T3_nmlk /= -2.0 * mat_P[m][l][k];

		return T3_nmlk;
	}

	ex fn_T4 (int n, int m, int l, int k){
		ex T4_nmlk;
  /**     Q_mlk + P_mlk.F_nmlk - beta_mlk.E_mlk)
		 ** T2_nmlk = ---------------------------------------- +
		 **  - 2.P_mlk
		 ** sqrt((Q_mlk+P_mlk.F_nmlk-beta_mlk.E_mlk)^2 - 4.P_mlk(Q_mlk.F_nmlk+beta_mlk.m_k^2-i.beta_mlk.rho))
		 **      - -------------------------------------------------------------------------------------------------
		 **- 2.P_mlk
   **/
		check0denom(mat_P[m][l][k], "T4", n, m, l, k);

		T4_nmlk = - ( mat_Q[m][l][k] + mat_P[m][l][k]*mat_F[n][m][l][k] - mat_beta[m][l][k]*mat_E[m][l][k]);
		T4_nmlk -= sqrt(
				pow(T4_nmlk,2)
				- 4.0 * mat_P[m][l][k] * (
						mat_Q[m][l][k]*mat_F[n][m][l][k]
						+ mat_beta[m][l][k]*mat_msquare[k]
						- I*mat_beta[m][l][k]*Rho /*This is Feynman's prescription*/
							 )
			       );
		T4_nmlk /= -2.0 * mat_P[m][l][k];
		return T4_nmlk;
	}
	
	void xloopsGiNaC_calc_lev4(){
		int k, l, m, n;
#ifdef _AGK_debugmode
		printf("Level 4 Calculating. Please Wait a Long Moment ...\n");
#endif
		// calc. z1phi, z2phi, ..., z1beta, z2beta, ...
		for(k=0; k<4; k++) for(l=0; l<4; l++) for(m = 0; m<4; m++){
			if(m!=l && m!=k && l!=k){
#ifdef _AGK_debugmode
				printf("%d-%d-%d\n", k+1, l+1, m+1);
#endif
				mat_z1phi[m][l][k] = fn_z1phi(m, l, k);
				mat_z2phi[m][l][k] = fn_z2phi(m, l, k);
				mat_z3phi[m][l][k] = fn_z3phi(m, l, k);
				mat_z4phi[m][l][k] = fn_z4phi(m, l, k);
				mat_z1beta[m][l][k] = fn_z1beta(m, l, k);
				mat_z2beta[m][l][k] = fn_z2beta(m, l, k);
				mat_z3beta[m][l][k] = fn_z3beta(m, l, k);
				mat_z4beta[m][l][k] = fn_z4beta(m, l, k);
			}
			else{ // set 0 value
				mat_z1phi[m][l][k] = 0;
				mat_z2phi[m][l][k] = 0;
				mat_z3phi[m][l][k] = 0;
				mat_z4phi[m][l][k] = 0;
				mat_z1beta[m][l][k] = 0;
				mat_z2beta[m][l][k] = 0;
				mat_z3beta[m][l][k] = 0;
				mat_z4beta[m][l][k] = 0;
			}
		}
		
		// calculate T1, T2, ...
		for(k=0; k<4; k++) for(l=0; l<4; l++) for(m=0; m<4; m++) for(n=0; n<4; n++){
			if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k){
#ifdef _AGK_debugmode
				printf("%d-%d-%d-%d\n", k+1, l+1, m+1, n+1);
#endif
				mat_T1[n][m][l][k] = fn_T1(n, m, l, k);
				mat_T2[n][m][l][k] = fn_T2(n, m, l, k);
				mat_T3[n][m][l][k] = fn_T3(n, m, l, k);
				mat_T4[n][m][l][k] = fn_T4(n, m, l, k);
			}
			else{ // set 0 value
				mat_T1[n][m][l][k] = 0;
				mat_T2[n][m][l][k] = 0;
				mat_T3[n][m][l][k] = 0;
				mat_T4[n][m][l][k] = 0;
			}
		}
	}
}// Namespace xloops
