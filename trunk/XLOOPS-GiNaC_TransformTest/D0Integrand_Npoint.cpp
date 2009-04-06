/*******************************************************************************
**
**      xloop-ginacs Project
**      For testing purpose, check NPoints equations
**
**
**      HCMUNS,
**      Mar 200
**
**
**      Author(s):      Son, D. H.      (sondo01@gmail.com)
**                      Khiem, Phan     (phanhongkhiem@gmail.com)
**                      Quan, Phan      (anhquan.phanle@gmail.com)
********************************************************************************
**
** Historial Log:
**      Date            Version Author          Description
**      _________       _______ _________       ________________________________
*******************************************************************************/

#include "ginac/ginac.h"

#include "trm.h"
#include "my_fns.h"
#include "trm2F.h"

#include "lev1.h"
#include "lev2.h"
#include "lev3.h"
#include "lev4.h"
#include "NPoint_Test.h"

using namespace GiNaC;

namespace xloops{

	/*
	To calculate equ 75, we need three things:
	1.	The constant
	2.	The positive integrand
	3.	The negative integrand
	*/
	ex NPoint_75_cte(int k, int l, int m, int n){
	/**
		The constant part of equ 75
	 */
		ex eq_75_cte = 0;

		return eq_75_cte;
	}

	ex NPoint_75_posIntegrand(int k, int l, int m, int n, const ex &z){
	/**
		The positive integrand of equ 75
	 */
		ex eq_75_posintegrand = 0;

		/* line 1*/
		eq_75_posintegrand = -(/* fg + (f-)g */ mat_f[l][k]*mat_g[m][l][k] + mat_fminus[l][k]*mat_g[m][l][k] )
			      * log(/*nom*/
				    (
				     - mat_P[m][l][k]/mat_beta[m][l][k] * pow(z,2)
				+ (mat_E[m][l][k] - mat_Q[m][l][k]/mat_beta[m][l][k]) * z
				- mat_msquare[k] + I*Rho
				    )
				/*denom*/
			/ (mat_Q[m][l][k] + mat_P[m][l][k]*z)
				   );

		/* line 2*/
		eq_75_posintegrand += mat_f[l][k]*mat_g[m][l][k]
				    * log(/*nom*/
					  (
					   - mat_P[m][l][k]*mat_phi[m][l][k] * pow(z,2)
				+ (mat_E[m][l][k] - mat_Q[m][l][k]*mat_phi[m][l][k]) * z
				- mat_msquare[k] + I*Rho
					  )
				/*denom*/
			/ (mat_Q[m][l][k] + mat_P[m][l][k]*z)
					 );

		/* line 3*/
		eq_75_posintegrand += mat_f[l][k]*mat_gminus[m][l][k]
					  * log(/*nom*/
				- (
				   - mat_P[m][l][k]*mat_phi[m][l][k] * pow(z,2)
				+ (mat_E[m][l][k] - mat_Q[m][l][k]*mat_phi[m][l][k]) * z
				- mat_msquare[k] + I*Rho
				  )
				/*denom*/
			/ (mat_Q[m][l][k] + mat_P[m][l][k]*z)
					       );
		/// mul with G(z)
		eq_75_posintegrand *= fn_G(k, l, m, n, z);

		return eq_75_posintegrand;
	}

	ex NPoint_75_negIntegrand(int k, int l, int m, int n, const ex &z){
	/**
		The negative integrand of equ 75
	 */
		ex eq_75_negintegrand = 0;


		/* line 1*/
		eq_75_negintegrand = mat_fminus[l][k]*mat_gminus[m][l][k]
				    * log(/*nom*/
					  (
					   - mat_P[m][l][k]/mat_beta[m][l][k] * pow(z,2)
				+ (mat_E[m][l][k] - mat_Q[m][l][k]/mat_beta[m][l][k]) * z
				- mat_msquare[k] + I*Rho
					  )
				/*denom*/
			/ (mat_Q[m][l][k] + mat_P[m][l][k]*z)
					 );
		/* line 2*/
		eq_75_negintegrand += mat_f[l][k]*mat_gminus[m][l][k]
					  * log(/*nom*/
				- (
				   - mat_P[m][l][k]/mat_beta[m][l][k] * pow(z,2)
				+ (mat_E[m][l][k] - mat_Q[m][l][k]/mat_beta[m][l][k]) * z
				- mat_msquare[k] + I*Rho
				  )
				/*denom*/
			/ (mat_Q[m][l][k] + mat_P[m][l][k]*z)
					       );

		/* line 3*/
		eq_75_negintegrand = -(/* (f-)(g-) + (f-)g */ (mat_fminus[l][k]*mat_gminus[m][l][k]) + mat_fminus[l][k]*mat_g[m][l][k])
					  * log(/*nom*/
				(
				 - mat_P[m][l][k]*mat_phi[m][l][k] * pow(z,2)
				+ (mat_E[m][l][k] - mat_Q[m][l][k]*mat_phi[m][l][k]) * z
				- mat_msquare[k] + I*Rho
				)
				/*denom*/
			/ (mat_Q[m][l][k] + mat_P[m][l][k]*z)
					       );

		/// mul with G(z)
		eq_75_negintegrand *= fn_G(k, l, m, n, z);

		return eq_75_negintegrand;
	}

	ex NPoint_75_01(){
/*
		return the positive integrand. May include the constant or not.
*/
		xloopsGiNaC_calc_lev1();
		xloopsGiNaC_calc_lev2();
		xloopsGiNaC_calc_lev3();
		xloopsGiNaC_calc_lev4();

		ex z = realsymbol("z"), t = realsymbol("t");
		ex integrand_75_01 = 0;

		int n, m, l, k;
		for (k=0; k<4; k++) for (l=0; l<4; l++) for (m=0; m<4; m++) for (n=0; n<4; n++){
			if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k){
				integrand_75_01 += NPoint_75_posIntegrand(k, l, m, n, z);
			}
		}

		return integrand_75_01;
	}

	ex NPoint_75_02(){
/*
		return the positive integrand. May include the constant or not.
*/
		xloopsGiNaC_calc_lev1();
		xloopsGiNaC_calc_lev2();
		xloopsGiNaC_calc_lev3();
		xloopsGiNaC_calc_lev4();

		ex z = realsymbol("z"), t = realsymbol("t");
		ex integrand_75_02 = 0;

		int n, m, l, k;
		for (k=0; k<4; k++) for (l=0; l<4; l++) for (m=0; m<4; m++) for (n=0; n<4; n++){
			if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k){
				integrand_75_02 += NPoint_75_negIntegrand(k, l, m, n, z);
			}
		}

		return integrand_75_02;
	}


ex NPoint_102_cte(int n, int m, int l, int k){
	ex factor; // the return factor
	const ex F = mat_F[n][m][l][k], beta = mat_beta[m][l][k], phi = mat_phi[m][l][k]
		T1 = fn_T1(n, m, l, k), T2 = fn_T2(n, m, l, k), T3 = fn_T3(n, m, l, k), T4 = fn_T4(n, m, l, k),
		flk = mat_f[l][k], gmlk = mat_g[m][l][k], fminuslk = mat_fminus[l][k], gminusmlk = mat_gminus[m][l][k];
	const ex Z1beta = fn_Z1(n, m, l, k, beta), Z2beta = fn_Z2(n, m, l, k, beta), Z1phi = fn_Z1(n, m, l, k, phi), Z2phi = fn_Z2(n, m, l, k, phi);

	ex factor1, factor2, factor3, factor4;

	//1
	factor1 = round_plus* abs(1.0 - mat_beta[m][l][k]*mat_phi[m][l][k])/ mat_P[m][l][k];
	//2
	factor2 = (flk+fminuslk) * gmlk * ( log(F/beta) - fn_Eta(n, m, l, k, beta) ) * fn_GZ(T3, T4);
	//3
	factor2+= flk*(gmlk + gminusmlk) * fn_Eta(n, m, l, k, phi) * fn_GZ(T3, T4);
	//4
	factor2+= gminusmlk*(
	                      (flk+fminuslk) * fn_Eta(n, m, l, k, beta) - fminuslk*log(F/beta) - flk*log(-F/beta)
	                    ) * fn_GZ(T1, T2);
	//5
	factor2+= -flk*(gmlk + gminusmlk) * fn_Eta(n, m, l, k, phi) * fn_GZ(T1, T2);
	//6
	factor2+= -flk*gmlk*fn_varLplus( (1.0 - beta*phi) / beta , F/beta, Rho /*use rho as the epsilon for now*/)
		- flk*gminusmlk*fn_varLplus( -(1.0 - beta*phi)/beta, -F/beta, Rho );
	//7
	factor2+= fminuslk*(gmlk + gminusmlk) * fn_varLminus((1.0 - beta*phi)/beta, F/beta, Rho);
	//8
	factor2+= -(flk+fminuslk)*gmlk * ( fn_varLplus(P*beta, -P*beta*Z1beta) + fn_varLplus(1.0, -Z2beta) );
	//9
	factor2+= flk*(gmlk + gminusmlk) * ( fn_varLplus(P, Q) + fn_varLminus(P, Q) )
		- fminuslk*(gmlk + gminusmlk) * ( fn_varLplus(P, Q) + fn_varLminus(P, Q) )

	return factor;
}

ex NPoint_102_posIntegrand(int n, int m, int l, int k, const ex &z){

}

ex NPoint_102_negIntegrand(int n, int m, int l, int k, const ex &z){
}

ex int102(){// integrand102
	ex return_int102 = 0;



	return return_int102;
} // integrand102

        /*******************************************************************************
	**      the main call from outside to calculate the integral
	********************************************************************************/
	ex NPoint_Equation(int eqnum, const ex& q_, const ex& m_, const ex& Rho_, const ex& Rho1_, const ex& Rho2_){
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
			case 7501: return_value = NPoint_75_01(); break;
			case 7502: return_value = NPoint_75_02(); break;
			default: return_value = 0;
		}

		return return_value;
	}

}// Namespace xloops

