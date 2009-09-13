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

	ex fn_z1beta (int m, int l, int k){
		ex z1beta_mlk;
	// init
//		ex E = fn_E(m, l, k), phi = fn_phi(m, l, k), P = fn_P(m, l, k), Q = fn_Q(m, l, k), msquare_k = mat_msquare[k], beta = fn_beta(m, l, k);
		ex E = evalf(mat_E[m][l][k]), phi = evalf(mat_phi[m][l][k]), P = evalf(mat_P[m][l][k]), Q = evalf(mat_Q[m][l][k]), msquare_k = mat_msquare[k], beta = evalf(mat_beta[m][l][k]);
		check0denom(beta, "z1beta", m, l, k);
		check0denom(P, "z1beta", m, l, k);

// calculation
		z1beta_mlk = - E + Q/beta;
		z1beta_mlk += sqrt(
				   pow(z1beta_mlk, 2) - 4.0*P*(msquare_k - I*Rho/*Feynman's pres.*/)/beta
				  );
		z1beta_mlk /= - 2.0 * P / beta;
		return z1beta_mlk;
	}
	
	ex fn_z2beta (int m, int l, int k){
		ex z2beta_mlk;
	// init
//		ex E = fn_E(m, l, k), phi = fn_phi(m, l, k), P = fn_P(m, l, k), Q = fn_Q(m, l, k), msquare_k = mat_msquare[k], beta = fn_beta(m, l, k);
		ex E = mat_E[m][l][k], phi = mat_phi[m][l][k], P = mat_P[m][l][k], Q = mat_Q[m][l][k], msquare_k = mat_msquare[k], beta = mat_beta[m][l][k];
		check0denom(beta, "z2beta", m, l, k);
		check0denom(P, "z2beta", m, l, k);

		// calculation
		z2beta_mlk = - E + Q/beta;
		z2beta_mlk -= sqrt(
				   pow(z2beta_mlk, 2) - 4.0*P*(msquare_k - I*Rho/*Feynman's pres.*/)/beta
				  );
		z2beta_mlk /= - 2.0 * P/beta;
		return z2beta_mlk;
	}

	ex fn_z1phi (int m, int l, int k){
		ex z1phi_mlk;
	// init
//		ex E = fn_E(m, l, k), phi = fn_phi(m, l, k), P = fn_P(m, l, k), Q = fn_Q(m, l, k), msquare_k = mat_msquare[k];
		ex E = mat_E[m][l][k], phi = mat_phi[m][l][k], P = mat_P[m][l][k], Q = mat_Q[m][l][k], msquare_k = mat_msquare[k];

		ex denom = -2.0 * P * phi;

		check0denom(denom, "z1phi", m, l, k);

	// calculation
		z1phi_mlk = - E + Q*phi;
		z1phi_mlk += sqrt( pow(z1phi_mlk, 2) - 4.0*P*phi*(msquare_k - I*Rho/*Feynman's pres.*/));
		z1phi_mlk /= denom;

		return z1phi_mlk;
	}

	ex fn_z2phi (int m, int l, int k){
		ex z2phi_mlk;
	// init
//		ex E = fn_E(m, l, k), phi = fn_phi(m, l, k), P = fn_P(m, l, k), Q = fn_Q(m, l, k), msquare_k = mat_msquare[k];
		ex E = mat_E[m][l][k], phi = mat_phi[m][l][k], P = mat_P[m][l][k], Q = mat_Q[m][l][k], msquare_k = mat_msquare[k];

		ex denom =  -2.0 * P * phi;

		check0denom(denom, "z2phi", m, l, k);

	// calculation
		z2phi_mlk =  - E + Q*phi;
		z2phi_mlk -= sqrt(
				  pow(z2phi_mlk, 2) - 4.0*P*phi*(msquare_k - I*Rho /*Feynman's pres.*/)
				 );
		z2phi_mlk /= denom;

		return z2phi_mlk;
	}

	ex fn_T1 (int n, int m, int l, int k){
		ex T1_nmlk;
		// init
//		ex P = fn_P(m, l, k), Q = fn_Q(m, l, k), F = fn_F(n, m, l, k), beta = fn_beta(m, l, k), E = fn_E(m, l, k), msquare_k = mat_msquare[k];
		ex P = mat_P[m][l][k], Q = mat_Q[m][l][k], F = mat_F[n][m][l][k], beta = mat_beta[m][l][k], E = mat_E[m][l][k], msquare_k = mat_msquare[k];
		check0denom(P, "T1", n, m, l, k);

		// calculate
		T1_nmlk = Q + P*F - beta*E;
		T1_nmlk += sqrt(
				pow(T1_nmlk,2) - 4.0*P*( Q*F + beta*msquare_k - I*beta*Rho /*Feynman's pres.*/)
			       );
		T1_nmlk /= -2.0 * P;

		return T1_nmlk;
	}

	ex fn_T2 (int n, int m, int l, int k){
		ex T2_nmlk;
		// init
//		ex P = fn_P(m, l, k), Q = fn_Q(m, l, k), F = fn_F(n, m, l, k), beta = fn_beta(m, l, k), E = fn_E(m, l, k), msquare_k = mat_msquare[k];
		ex P = mat_P[m][l][k], Q = mat_Q[m][l][k], F = mat_F[n][m][l][k], beta = mat_beta[m][l][k], E = mat_E[m][l][k], msquare_k = mat_msquare[k];
		check0denom(P, "T2", n, m, l, k);

		T2_nmlk = Q + P*F - beta*E;
		T2_nmlk -= sqrt(
				pow(T2_nmlk,2) - 4.0*P*( Q*F + beta*msquare_k - I*beta*Rho /*Feynman's pres.*/)
			       );
		T2_nmlk /= -2.0*P;
		return T2_nmlk;
	}


	void lev4Calc(){
		// calculate previous level
		lev3Calc();
		int n, m, l, k;
		for(k=0; k<4; k++) for(l=0; l<4; l++) for(m=0; m<4; m++) if(m!=l && m!=k && l!=k){
			mat_z1phi[m][l][k] = fn_beta(m, l, k);
			mat_z2phi[m][l][k] = fn_phi(m, l, k);
			mat_z1beta[m][l][k] = fn_g(m, l, k);
			mat_z2beta[m][l][k] = fn_gminus(m, l, k);
		}
		for(k=0; k<4; k++) for(l=0; l<4; l++) for(m=0; m<4; m++) for(n=0; n<4; n++) if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k){
			mat_T1[n][m][l][k] = fn_T1(n, m, l, k);
			mat_T2[n][m][l][k] = fn_T2(n, m, l, k);
		}
	}
	void lev4Calc(int debug){
		if (debug == 0 ){
			lev4Calc();
			return;
		}
		// calculate previous level
		lev3Calc(1);
		printf("Level 3 calculations...\n");

		int n, m, l, k;
		for(k=0; k<4; k++) for(l=0; l<4; l++) for(m=0; m<4; m++) if(m!=l && m!=k && l!=k){
			printf("z1phi[%d,%d,%d]\n", m, l, k);
			mat_z1phi[m][l][k] = fn_beta(m, l, k);
			printf("z2phi[%d,%d,%d]\n", m, l, k);
			mat_z2phi[m][l][k] = fn_phi(m, l, k);
			printf("z1beta[%d,%d,%d]\n", m, l, k);
			mat_z1beta[m][l][k] = fn_g(m, l, k);
			printf("z2beta[%d,%d,%d]\n", m, l, k);
			mat_z2beta[m][l][k] = fn_gminus(m, l, k);
		}
		for(k=0; k<4; k++) for(l=0; l<4; l++) for(m=0; m<4; m++) for(n=0; n<4; n++) if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k){
			printf("T1[%d,%d,%d,%d]\n", n, m, l, k);
			mat_T1[n][m][l][k] = fn_T1(n, m, l, k);
			printf("T2[%d,%d,%d,%d]\n", n, m, l, k);
			mat_T2[n][m][l][k] = fn_T2(n, m, l, k);
		}
	}
}// Namespace xloops
