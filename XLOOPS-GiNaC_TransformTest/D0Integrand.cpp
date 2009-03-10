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
**	Calculate the Integrands
**
********************************************************************************
**
** Historial Log:
**	Date		Version	Author		Description
**	_________	_______	_________	________________________________
*******************************************************************************/
//#define _AGK_debugmode

#include "ginac/ginac.h"

#include "trm.h"
#include "my_fns.h"
#include "trm2F.h"

#include "lev1.h"
#include "lev2.h"
#include "lev3.h"
#include "lev4.h"
#include "D0Integrand.h"

using namespace GiNaC;

namespace xloops{

ex D0_integrand1(){
	xloopsGiNaC_calc_lev1();
	// integrand1 = 1/(P1*P2*P3*P4)
	ex l0 = realsymbol("l0"), l1 = realsymbol("l1"), l2 = realsymbol("l2"), lortho = realsymbol("lortho");

	ex P1, P2, P3, P4;
	ex integrand1;

	P1 = pow(l0 + mat_q[0][0], 2) - pow(l1, 2) - pow(l2, 2) - pow(lortho, 2) - mat_msquare[0] + I*Rho;
	P2 = pow(l0 + mat_q[1][0], 2) - pow(l1 + mat_q[1][1], 2) - pow(l2, 2) - pow(lortho, 2) - mat_msquare[1] + I*Rho;
	P3 = pow(l0 + mat_q[2][0], 2) - pow(l1 + mat_q[2][1], 2) - pow(l2 + mat_q[2][2], 2) - pow(lortho, 2) - mat_msquare[2] + I*Rho;
	P4 = pow(l0, 2) - pow(l1, 2) - pow(l2, 2) - pow(lortho, 2) - mat_msquare[3] + I*Rho;

	integrand1 = 1.0 / (P1*P2*P3*P4);

	return integrand1;
}

ex D0_integrand9(){
	xloopsGiNaC_calc_lev1();

		// integrand9 = 1/(l0^2 - l1^2 -...) / PI(alk.l0 + blk.l1 + ...)
	ex integrand9 = 0;

	ex l0 = realsymbol("l0"), l1 = realsymbol("l1"), l2 = realsymbol("l2"), lortho = realsymbol("lortho");

	int l, k;

	for(k = 0; k<4; k++){
		ex term1 = (pow(l0, 2) - pow(l1, 2) - pow(l2, 2) - pow(lortho, 2) - mat_msquare[k] + I*Rho);
		ex term2 = 1;
		for(l = 0; l<4; l++){
			if(l!=k){
				term2 *=(mat_a[l][k]*l0 + mat_b[l][k]*l1 + mat_c[l][k]*l2 + mat_d[l][k]);
			}
		}
		integrand9 += (1.0/term1) * (1.0/term2);
	}

	return integrand9;
}

ex D0_integrand12(){
	xloopsGiNaC_calc_lev1();
	xloopsGiNaC_calc_lev2(); // AC is in level 2 terms

	ex integrand12 = 0;

	ex x = realsymbol("x"), y = realsymbol("y"), z = realsymbol("z"), t = realsymbol("t");

	int l, k;

	for(k = 0; k<4; k++){
		ex term1 = 2.0*x*z - pow(z, 2) - pow(y, 2) - pow(t, 2) - mat_msquare[k] + I*Rho;

		ex term2 = 1;
		for(l = 0; l<4; l++){
			if(l!=k){
				term2 *=(mat_a[l][k]*z + mat_b[l][k]*y + mat_AC[l][k]*x + mat_d[l][k]);
			}
		}

		integrand12 += (1.0/term1) * (1.0/term2);
	}

	return integrand12;
}

ex D0_integrand1801(){ // equation 18, D0+ (01)
	xloopsGiNaC_calc_lev1();
	xloopsGiNaC_calc_lev2(); // AC_lk in lev 2 calculation

	ex integrand1801 = 0;

	ex x = realsymbol("x"), y = realsymbol("y"), z = realsymbol("z"), t = realsymbol("t");

	int m, l, k;

	for(k = 0; k<4; k++) for(l = 0; l<4; l++){
		if(l!=k){
			ex term1 = 1.0/mat_AC[l][k];
			ex term2 = 1;
			for(m=0; m<4; m++){
				if(m!=l && m!=k){
					term2 *= 1.0/(mat_A[m][l][k]*z + mat_B[m][l][k]*y + mat_C[m][l][k]);
					term2 *= mat_f[l][k] * (1.0 - myfn_delta(mat_AC[l][k]));
					term2 *= 1.0/( (1.0-2.0*mat_a[l][k]/mat_AC[l][k])*pow(z,2)
					               - 2.0*mat_b[l][k]*y*z/mat_AC[l][k]
					               - 2.0*mat_d[l][k]*z/mat_AC[l][k]
					               -pow(y, 2)
					               -pow(t, 2)
					               -mat_msquare[k]
					               + I*Rho);
				}
			}
			integrand1801 += term1 * term2;
		}
	}

		// because we get the integral for t limit = -infty..infty
		//	we will drop the *2 factor
	integrand1801 *= Pi*I;
	return integrand1801;
}

ex D0_integrand1802(){ // equation 18, D0-
	xloopsGiNaC_calc_lev1();
	xloopsGiNaC_calc_lev2(); // AC_lk in lev 2 calculation

	ex integrand1802 = 0;

	ex x = realsymbol("x"), y = realsymbol("y"), z = realsymbol("z"), t = realsymbol("t");

	int m, l, k;

	for(k = 0; k<4; k++) for(l = 0; l<4; l++){
		if(l!=k){
			ex term1 = 1.0/mat_AC[l][k];
			ex term2 = 1;
			for(m=0; m<4; m++){
				if(m!=l && m!=k){
					term2 *= 1.0/(mat_A[m][l][k]*z + mat_B[m][l][k]*y + mat_C[m][l][k]);
						/*the only different between D0+ and - is mat_f and mat_fminus*/
					term2 *= mat_fminus[l][k] * (1.0 - myfn_delta(mat_AC[l][k]));
					term2 *= 1.0/( (1.0-2.0*mat_a[l][k]/mat_AC[l][k])*pow(z,2)
					               - 2.0*mat_b[l][k]*y*z/mat_AC[l][k]
					               - 2.0*mat_d[l][k]*z/mat_AC[l][k]
					               -pow(y, 2)
					               -pow(t, 2)
					               -mat_msquare[k]
					               + I*Rho);
				}
			}
			integrand1802 += term1 * term2;
		}
	}

		// the integrate is for limit of t = -infty..infty => no *2 factor
	integrand1802 *= - Pi*I;
	return integrand1802;
}


ex D0_integrand2501(){ // equation 25, D0+ (01)
	xloopsGiNaC_calc_lev1();
	xloopsGiNaC_calc_lev2(); // AC_lk in lev 2 calculation

	ex integrand2501 = 0;

	ex x = realsymbol("x"), y = realsymbol("y"), z = realsymbol("z"), t = realsymbol("t");

	int m, l, k;

	for(k = 0; k<4; k++) for(l = 0; l<4; l++){
		if(l!=k){
			ex term1 = 1.0/mat_AC[l][k];
			ex term2 = 1;
			for(m=0; m<4; m++){
				if(m!=l && m!=k){
					term2 *= 1.0/(mat_A[m][l][k]*z + mat_B[m][l][k]*y + mat_C[m][l][k]);
					term2 *= mat_f[l][k] * (1.0 - myfn_delta(mat_AC[l][k]));
					term2 *= 1.0/( (1.0-2.0*mat_a[l][k]/mat_AC[l][k])*pow(z,2)
					               - 2.0*mat_b[l][k]*y*z/mat_AC[l][k]
					               - 2.0*mat_d[l][k]*z/mat_AC[l][k]
					               -pow(y, 2)
					               +pow(t, 2)
					               -mat_msquare[k]
					               + I*Rho);
				}
			}
			integrand2501 += term1 * term2;
		}
	}

	integrand2501 *= Pi;
	return integrand2501;
}

ex D0_integrand2502(){ // equation 25, D0-
	xloopsGiNaC_calc_lev1();
	xloopsGiNaC_calc_lev2(); // AC_lk in lev 2 calculation

	ex integrand2502 = 0;

	ex x = realsymbol("x"), y = realsymbol("y"), z = realsymbol("z"), t = realsymbol("t");

	int m, l, k;

	for(k = 0; k<4; k++) for(l = 0; l<4; l++){
		if(l!=k){
			ex term1 = 1.0/mat_AC[l][k];
			ex term2 = 1;
			for(m=0; m<4; m++){
				if(m!=l && m!=k){
					term2 *= 1.0/(mat_A[m][l][k]*z + mat_B[m][l][k]*y + mat_C[m][l][k]);
						/*the only different between D0+ and - is mat_f and mat_fminus*/
					term2 *= mat_fminus[l][k] * (1.0 - myfn_delta(mat_AC[l][k]));
					term2 *= 1.0/( (1.0-2.0*mat_a[l][k]/mat_AC[l][k])*pow(z,2)
					               - 2.0*mat_b[l][k]*y*z/mat_AC[l][k]
					               - 2.0*mat_d[l][k]*z/mat_AC[l][k]
					               -pow(y, 2)
					               +pow(t, 2)
					               -mat_msquare[k]
					               + I*Rho);
				}
			}
			integrand2502 += term1 * term2;
		}
	}

		// the integrate is for limit of t = -infty..infty => no *2 factor
	integrand2502 *= - Pi;
	return integrand2502;
}

	/** For integral 30, we need to calculate I_nmlk
	 */
ex fn_I(int k, int l, int m, int n, const ex &z, const ex &t){
	ex I_nmlk;

	I_nmlk = 1.0/mat_AC[l][k];
	I_nmlk *= (1.0 - myfn_delta(mat_AC[l][k]))*(1.0 - myfn_delta(mat_B[m][l][k])) / (mat_A[n][l][k]*mat_B[m][l][k] - mat_A[m][l][k]*mat_B[n][l][k]);
	I_nmlk *= 1.0/(z + mat_F[n][m][l][k]);
	I_nmlk *= 1.0/( mat_D[m][l][k]*pow(z,2)
	                - 2.0*(mat_A[m][l][k]*mat_B[m][l][k] - mat_alpha[l][k])*z*t
	                - 2.0*mat_C[m][l][k]*t/mat_B[m][l][k]
	                - 2.0*mat_d[l][k]*z/mat_AC[l][k]
	                + pow(t,2)
	                - mat_msquare[k]
	                + I*Rho
	              );

	return I_nmlk;
}
	/** Now we can go calculating the integral 30 (there are 4 sub-integral), and I call them Integrand3001, 3002, 3003, and 3004.
	 */
ex D0_integrand3001(){
	xloopsGiNaC_calc_lev1();
	xloopsGiNaC_calc_lev2();
	xloopsGiNaC_calc_lev3();

	ex z = realsymbol("z"), t = realsymbol("t");
	ex integrand3001 = 0;

	int n, m, l, k;
	for (k=0; k<4; k++) for (l=0; l<4; l++) for (m=0; m<4; m++) for (n=0; n<4; n++){
		if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k){
			integrand3001 += mat_f[l][k]*mat_g[m][l][k]*fn_I(n, m, l, k, z, t);
		}
	}

	integrand3001 *= I*pow(Pi,2);
	return integrand3001;
}

ex D0_integrand3002(){
	xloopsGiNaC_calc_lev1();
	xloopsGiNaC_calc_lev2();
	xloopsGiNaC_calc_lev3();

	ex z = realsymbol("z"), t = realsymbol("t");
	ex integrand3002 = 0;

	int n, m, l, k;
	for (k=0; k<4; k++) for (l=0; l<4; l++) for (m=0; m<4; m++) for (n=0; n<4; n++){
		if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k){
			integrand3002 += mat_f[l][k]*mat_gminus[m][l][k]*fn_I(n, m, l, k, z, t);
		}
	}

	integrand3002 *= -I*pow(Pi,2);
	return integrand3002;
}

ex D0_integrand3003(){
	xloopsGiNaC_calc_lev1();
	xloopsGiNaC_calc_lev2();
	xloopsGiNaC_calc_lev3();

	ex z = realsymbol("z"), t = realsymbol("t");
	ex integrand3003 = 0;

	int n, m, l, k;
	for (k=0; k<4; k++) for (l=0; l<4; l++) for (m=0; m<4; m++) for (n=0; n<4; n++){
		if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k){
			integrand3003 += mat_fminus[l][k]*mat_g[m][l][k]*fn_I(n, m, l, k, z, t);
		}
	}

	integrand3003 *= -I*pow(Pi,2);
	return integrand3003;
}

ex D0_integrand3004(){
	xloopsGiNaC_calc_lev1();
	xloopsGiNaC_calc_lev2();
	xloopsGiNaC_calc_lev3();

	ex z = realsymbol("z"), t = realsymbol("t");
	ex integrand3004 = 0;

	int n, m, l, k;
	for (k=0; k<4; k++) for (l=0; l<4; l++) for (m=0; m<4; m++) for (n=0; n<4; n++){
		if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k){
			integrand3004 += mat_fminus[l][k]*mat_gminus[m][l][k]*fn_I(n, m, l, k, z, t);
		}
	}

	integrand3004 *= -I*pow(Pi,2);
	return integrand3004;
}


	// to calculate integrand 40, we will need the function G(z)
ex fn_G(int k, int l, int m, int n, const ex &z){
	ex return_G = 0;

	return_G = mat_beta[m][l][k] * ( mat_E[m][l][k]*z - mat_msquare[k] + I*Rho );
	return_G += -( mat_P[m][l][k]*z + mat_Q[m][l][k]) * (z + mat_F[n][m][l][k]);
	return_G = 1.0/return_G;

	return return_G;
}

ex D0_integrand4001(){ // again, the integrand 40 has two parts z>=0 and z<0, 4001 is for z>=0
	xloopsGiNaC_calc_lev1();
	xloopsGiNaC_calc_lev2();
	xloopsGiNaC_calc_lev3();

	ex z = realsymbol("z");

	ex integrand4001 = 0;

	int n, m, l, k;
	for (k=0; k<4; k++) for (l=0; l<4; l++) for (m=0; m<4; m++) for (n=0; n<4; n++){
		if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k){
			// 1.	calculate the jacobi
			ex jacobi = 1.0/mat_AC[l][k];
			jacobi *= 1.0/(mat_B[m][l][k]*mat_A[n][l][k] - mat_B[n][l][k]*mat_A[m][l][k]);
			jacobi *= ( 1.0 - myfn_delta(mat_AC[l][k]) ) * ( 1.0 - myfn_delta(mat_B[m][l][k]) ) * abs( 1.0 - mat_beta[m][l][k]*mat_phi[m][l][k] );


			ex tmp_integrand = 0;
			// line 1
			tmp_integrand = ( /*  fg + (f-)g  */
			                  mat_f[l][k]*mat_g[m][l][k] + mat_fminus[l][k]*mat_g[m][l][k])
				* log(/*  F/beta  */
				      mat_F[n][m][l][k]/mat_beta[m][l][k]
				     );

			// line 2
			tmp_integrand += -mat_f[l][k]*mat_g[m][l][k]
				* log( (( 1.0 - mat_beta[m][l][k]*mat_phi[m][l][k] )*z + mat_F[n][m][l][k])
				       /
				       mat_beta[m][l][k])

				- mat_f[l][k]*mat_gminus[m][l][k]
				* log(-(( 1.0 - mat_beta[m][l][k]*mat_phi[m][l][k])*z + mat_F[n][m][l][k])
				      /
				      mat_beta[m][l][k]);

			// line 3
			tmp_integrand += - ( /*  fg + (f-)g  */
			                     mat_f[l][k]*mat_g[m][l][k] + mat_fminus[l][k]*mat_g[m][l][k])
				* log ( // begin log
					/* nominator */
				        ( -	(mat_P[m][l][k]/mat_beta[m][l][k]) * pow(z,2)
				          +	(mat_E[m][l][k] - mat_Q[m][l][k]/mat_beta[m][l][k]) * z
				          -	mat_msquare[k] + I*Rho)
					/* denominator */
				        /(mat_Q[m][l][k] + mat_P[m][l][k]*z)
				      ); // end log
			// line 4
			tmp_integrand += mat_f[l][k]*mat_g[m][l][k]
				* log ( // begin log
					/* nominator */
				        ( -	(mat_P[m][l][k]*mat_phi[m][l][k]) * pow(z,2)
				          +	(mat_E[m][l][k] - mat_Q[m][l][k]*mat_phi[m][l][k]) * z
				          -	mat_msquare[k] + I*Rho)
					/* denominator */
				        /(mat_Q[m][l][k] + mat_P[m][l][k]*z)
				      ); // end log
			// line 5
			tmp_integrand += mat_f[l][k]*mat_g[m][l][k]
				* log ( // begin log
					/* nominator */
				        (
				          -( -	(mat_P[m][l][k]*mat_phi[m][l][k]) * pow(z,2)
				             +	(mat_E[m][l][k] - mat_Q[m][l][k]*mat_phi[m][l][k]) * z
				             -	mat_msquare[k] + I*Rho
				           )
				        )
					/* denominator */
				        /(mat_Q[m][l][k] + mat_P[m][l][k]*z)
				      ); // end log

			// sum parts of the loop
			integrand4001 += jacobi * fn_G(k, l, m, n, z) * tmp_integrand;
		} // indices are different
	} // main for loop through indices

	integrand4001 *= I*pow(Pi,2);
	return integrand4001;
}

ex D0_integrand4002(){ // again, the integrand 40 has two parts z>=0 and z<0, 4002 is for z>=0
	xloopsGiNaC_calc_lev1();
	xloopsGiNaC_calc_lev2();
	xloopsGiNaC_calc_lev3();

	ex z = realsymbol("z");

	ex integrand4002;
	int n, m, l, k;
	for (k=0; k<4; k++) for (l=0; l<4; l++) for (m=0; m<4; m++) for (n=0; n<4; n++){
		if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k){
			// 1.	calculate the jacobi
			ex jacobi = 1.0/mat_AC[l][k];
			jacobi *= 1.0/(mat_B[m][l][k]*mat_A[n][l][k] - mat_B[n][l][k]*mat_A[m][l][k]);
			jacobi *= ( 1.0 - myfn_delta(mat_AC[l][k]) ) * ( 1.0 - myfn_delta(mat_B[m][l][k]) ) * abs( 1.0 - mat_beta[m][l][k]*mat_phi[m][l][k] );


			ex tmp_integrand = 0;
			// line 1
			tmp_integrand = -mat_fminus[l][k]*mat_gminus[m][l][k]
				* log( mat_F[n][m][l][k]/mat_beta[m][l][k] )
				- mat_f[l][k]*mat_gminus[m][l][k]
				* log( -mat_F[n][m][l][k]/mat_beta[m][l][k] );

			// line 2
			tmp_integrand += ( /*  (f-)(g-) + (f-)g  */
			                   mat_fminus[l][k]*mat_gminus[m][l][k] + mat_fminus[l][k]*mat_g[m][l][k]
			                 )
				* log( (( 1.0 - mat_beta[m][l][k]*mat_phi[m][l][k] )*z + mat_F[n][m][l][k])
				       /
				       mat_beta[m][l][k]);

			// line 3
			tmp_integrand += mat_fminus[l][k]*mat_gminus[m][l][k]
				* log ( // begin log
					/* nominator */
				        ( -	(mat_P[m][l][k]/mat_beta[m][l][k]) * pow(z,2)
				          +	(mat_E[m][l][k] - mat_Q[m][l][k]/mat_beta[m][l][k]) * z
				          -	mat_msquare[k] + I*Rho)
					/* denominator */
				        /(mat_Q[m][l][k] + mat_P[m][l][k]*z)
				      ); // end log
			// line 4
			tmp_integrand += mat_f[l][k]*mat_gminus[m][l][k]
				* log ( // begin log
					/* nominator */
				        (
				          -(-	(mat_P[m][l][k]/mat_beta[m][l][k]) * pow(z,2)
				            +	(mat_E[m][l][k] - mat_Q[m][l][k]/mat_beta[m][l][k]) * z
				            -	mat_msquare[k] + I*Rho
				           )
				        )
					/* denominator */
				        /(mat_Q[m][l][k] + mat_P[m][l][k]*z)
				      ); // end log
			// line 5
			tmp_integrand += - /*  (f-)(g-) + (f-)g  */
			( mat_fminus[l][k]*mat_gminus[m][l][k] + mat_fminus[l][k]*mat_g[m][l][k] )
				* log ( // begin log
					/* nominator */
				        (-	(mat_P[m][l][k]*mat_phi[m][l][k]) * pow(z,2)
				         +	(mat_E[m][l][k] - mat_Q[m][l][k]*mat_phi[m][l][k]) * z
				         -	mat_msquare[k] + I*Rho
				        )
					/* denominator */
				        /(mat_Q[m][l][k] + mat_P[m][l][k]*z)
				      ); // end log
			// sum parts of the loop
			integrand4002 += jacobi * fn_G(k, l, m, n, z) * tmp_integrand;
		} // indices are different
	} // main for loop through indices

	integrand4002 *= I*pow(Pi,2);
	return integrand4002;
}

//	For equation 48, we will need Coff calculate
ex D0_integrand48_Coff(int k, int l, int m, int n){
	ex Coff;

	// I split the Coff into 2 line, as in the Oneloop4pt.pdf
	Coff = (1.0/mat_AC[l][k]) * (  1.0/(mat_B[m][l][k]*mat_A[n][l][k] - mat_B[n][l][k]*mat_A[m][l][k])  );

	Coff *= (1.0-myfn_delta(mat_AC[l][k])) * (1.0-myfn_delta(mat_B[m][l][k])) * abs(1.0-mat_beta[m][l][k]*mat_phi[m][l][k]);

	return Coff;
}
//	for equation 48, we also need pos term, neg term and extra term
ex D0_integrand48_posTerm(int k, int l, int m, int n, const ex &z){
	ex posTerm;

	// line 1
	posTerm = mat_OPlus[n][m][l][k]
		- mat_f[l][k]*mat_g[m][l][k]*log(
		                                  ( (1.0-mat_beta[m][l][k]*mat_phi[m][l][k])/mat_beta[m][l][k] )*z
		                                  + mat_F[n][m][l][k]/mat_beta[m][l][k]
		                                );
	// line 2
	posTerm +=
		- mat_f[l][k]*mat_gminus[m][l][k]*log(
		                                       ( -(1.0-mat_beta[m][l][k]*mat_phi[m][l][k])/mat_beta[m][l][k] )*z
		                                       - mat_F[n][m][l][k]/mat_beta[m][l][k]
		                                     );
	// line 3
	posTerm += /*fg + (f-)g*/
	- (mat_f[l][k]*mat_g[m][l][k]+mat_fminus[l][k]*mat_g[m][l][k]) *log(
		( -mat_P[m][l][k]*z / mat_beta[m][l][k] + mat_P[m][l][k]*mat_z1beta[m][l][k] / mat_beta[m][l][k])
		);
	// line 4
	posTerm += /*fg + (f-)g*/ - (mat_f[l][k]*mat_g[m][l][k]+mat_fminus[l][k]*mat_g[m][l][k]) *log( z-mat_z2beta[m][l][k] );
	// line 5
	posTerm += mat_f[l][k]*mat_g[m][l][k] * log(-mat_P[m][l][k]*mat_phi[m][l][k]*z + mat_P[m][l][k]*mat_phi[m][l][k]*mat_z1phi[m][l][k]) + mat_f[l][k]*mat_g[m][l][k]*log(z + mat_z2phi[m][l][k]);
	// line 6
	posTerm += mat_f[l][k]*mat_gminus[m][l][k] * log(mat_P[m][l][k]*mat_phi[m][l][k]*z - mat_P[m][l][k]*mat_phi[m][l][k]*mat_z1phi[m][l][k]) + mat_f[l][k]*mat_g[m][l][k]*log(z + mat_z2phi[m][l][k]);
	// line 7
	posTerm += /*f(g-) + (f-)g*/ (mat_fminus[l][k]*mat_g[m][l][k]+mat_f[l][k]*mat_gminus[m][l][k]) *log( mat_P[m][l][k]*z + mat_Q[m][l][k] );

	posTerm *= fn_G(k, l, m, n, z);
	return posTerm;
}

//	for equation 48, we also need pos term, neg term and extra term
ex D0_integrand48_negTerm(int k, int l, int m, int n, const ex &z){
	ex negTerm;

	// line 1
	negTerm = mat_OMinus[n][m][l][k]
		+ (mat_fminus[l][k]*mat_gminus[m][l][k] + mat_fminus[l][k]*mat_g[m][l][k])*log(
		                                  ( (1.0-mat_beta[m][l][k]*mat_phi[m][l][k])/mat_beta[m][l][k] )*z
		                                  + mat_F[n][m][l][k]/mat_beta[m][l][k]
		                                );
	// line 2
	negTerm +=
		mat_fminus[l][k]*mat_gminus[m][l][k]*log(
		                                          - (-mat_P[m][l][k]*z/mat_beta[m][l][k])
		                                          + mat_P[m][l][k]*mat_z1beta[m][l][k]/mat_beta[m][l][k]
		                                     )
		+	mat_fminus[l][k]*mat_gminus[m][l][k]*log(z-mat_z2beta[m][l][k]);
	// line 3
	negTerm += mat_f[l][k]*mat_gminus[m][l][k]*log(
	                                                -mat_P[m][l][k]*z/mat_beta[m][l][k]
	                                                -mat_P[m][l][k]*mat_z1beta[m][l][k]/mat_beta[m][l][k]
	                                              )
		+ mat_fminus[l][k]*mat_gminus[m][l][k] * log (z-mat_z2beta[m][l][k]);
	// line 4
	negTerm += - (mat_fminus[l][k]*mat_gminus[m][l][k]+mat_fminus[l][k]*mat_g[m][l][k]) *log( -mat_P[m][l][k]*mat_phi[m][l][k]*z + mat_P[m][l][k]*mat_phi[m][l][k]*mat_z1phi[m][l][k] )
		- (mat_fminus[l][k]*mat_gminus[m][l][k] + mat_fminus[l][k]*mat_g[m][l][k]) * log (z - mat_z2phi[m][l][k]);
	// line 5
	negTerm += (mat_fminus[l][k]*mat_gminus[m][l][k] - mat_f[l][k]*mat_gminus[m][l][k]) * log (mat_P[m][l][k]*z + mat_Q[m][l][k]);

	negTerm *= fn_G(k, l, m, n, z);
	return negTerm;
}

ex D0_integrand48_extraTerm(int k, int l, int m, int n, const ex &z){
	ex extraTerm;

	ex A0 = imag_part(mat_P[m][l][k]*mat_E[m][l][k]);
	ex B0 = imag_part(mat_E[m][l][k] - mat_P[m][l][k]*mat_msquare[k] + I*Rho*mat_P[m][l][k]);
	ex C0 = imag_part((-mat_msquare[k] + I*Rho)*mat_Q[k][l][m].conjugate());

	extraTerm = 2.0*Pi*I*mat_fminus[l][k]*mat_g[m][l][k] * my_step(-imag_part(mat_Q[m][l][k])) * my_step(A0*pow(z,2) + B0*z + C0) * fn_G(k, l, m, n, z);
	extraTerm += 2.0*Pi*I*mat_f[l][k]*mat_gminus[m][l][k] * my_step(imag_part(mat_Q[m][l][k])) * my_step(-A0*pow(z,2) - B0*z - C0) * fn_G(k, l, m, n, z);

	return extraTerm;
}

ex D0_integrand4801(){
	xloopsGiNaC_calc_lev1();
	xloopsGiNaC_calc_lev2();
	xloopsGiNaC_calc_lev3();

	ex z = realsymbol("z");
	ex integrand4801 = 0;
	int n, m, l, k;
	for (k=0; k<4; k++) for (l=0; l<4; l++) for (m=0; m<4; m++) for (n=0; n<4; n++){
		if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k){
			integrand4801 += D0_integrand48_Coff(k, l, m ,n) * D0_integrand48_posTerm(k, l, m, n, z);
		}
	}

	integrand4801 *= I*pow(Pi, 2);

	return integrand4801;
}


ex D0_integrand4802(){
	xloopsGiNaC_calc_lev1();
	xloopsGiNaC_calc_lev2();
	xloopsGiNaC_calc_lev3();

	ex z = realsymbol("z");
	ex integrand4802 = 0;
	int n, m, l, k;
	for (k=0; k<4; k++) for (l=0; l<4; l++) for (m=0; m<4; m++) for (n=0; n<4; n++){
		if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k){
			integrand4802 += D0_integrand48_Coff(k, l, m ,n) * D0_integrand48_negTerm(k, l, m, n, z);
		}
	}

	integrand4802 *= I*pow(Pi, 2);

	return integrand4802;
}

ex D0_integrand4803(){
	xloopsGiNaC_calc_lev1();
	xloopsGiNaC_calc_lev2();
	xloopsGiNaC_calc_lev3();

	ex z = realsymbol("z");
	ex integrand4803 = 0;
	int n, m, l, k;
	for (k=0; k<4; k++) for (l=0; l<4; l++) for (m=0; m<4; m++) for (n=0; n<4; n++){
		if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k){
			integrand4803 += D0_integrand48_extraTerm(k, l, m, n, z);
		}
	}

	integrand4803 *= I*pow(Pi, 2);

	return integrand4803;
}

	/*******************************************************************************
	**	the main call from outside to calculate the integral
	********************************************************************************/
ex D0_integrand(int eqnum, const ex& q_, const ex& m_, const ex& Rho_, const ex& Rho1_, const ex& Rho2_){
	ex return_value;

		// init variables
	mat_msquare[0]=m_.op(0); mat_msquare[1]=m_.op(1); mat_msquare[2]=m_.op(2); mat_msquare[3]=m_.op(3);
	mat_q[0][0]=q_.op(0);mat_q[0][1]=0;mat_q[0][2]=0;mat_q[0][3]=0;
	mat_q[1][0]=q_.op(1);mat_q[1][1]=q_.op(2);mat_q[1][2]=0;mat_q[1][3]=0;
	mat_q[2][0]=q_.op(3);mat_q[2][1]=q_.op(4);mat_q[2][2]=q_.op(5);mat_q[2][3]=0;
	input_Rho = Rho_; input_Rho1 = Rho1_; input_Rho2 = Rho2_;
	Rho = Rho_; Rho1 = Rho1_; Rho2 = Rho2_;

		// return the desire integrand;
	switch(eqnum){
	case 1: return_value = D0_integrand1(); break;
	case 9: return_value = D0_integrand9(); break;
	case 12: return_value = D0_integrand12(); break;
	case 1801/*D0+*/: return_value = D0_integrand1801(); break;
	case 1802/*D0-*/: return_value = D0_integrand1802(); break;
	case 3001/*D0++*/: return_value = D0_integrand3001(); break;
	case 3002/*D0+-*/: return_value = D0_integrand3002(); break;
	case 3003/*D0-+*/: return_value = D0_integrand3003(); break;
	case 3004/*D0--*/: return_value = D0_integrand3004(); break;
	case 4001/*D0+*/: return_value = D0_integrand4001(); break;
	case 4002/*D0-*/: return_value = D0_integrand4002(); break;
	case 4801/*D0-*/: return_value = D0_integrand4801(); break;
	case 4802/*D0-*/: return_value = D0_integrand4802(); break;
	case 4803/*D0-*/: return_value = D0_integrand4803(); break;
	default: return_value = 0;
	}

	return return_value;
}
}// Namespace xloops
