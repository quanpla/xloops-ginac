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
	/**
		G(z) & GZ(z)
	*/
	ex fn_G(int k, int l, int m, int n, const ex &z){
		ex return_G = 0;

		return_G = mat_beta[m][l][k] * ( mat_E[m][l][k]*z - mat_msquare[k] + I*Rho );
		return_G += -( mat_P[m][l][k]*z + mat_Q[m][l][k]) * (z + mat_F[n][m][l][k]);
		return_G = 1.0/return_G;

		return return_G;
	}
	
	ex fn_GZ(const ex &x, const ex &y){
		ex return_GZ = 0;
		
		return_GZ = - (log(x) - log(y))/(x - y);
		
		return return_GZ;
	}
	
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

