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

#include <iostream>

#include "ginac/ginac.h"

#include "extrm.h"
#include "my_fns.h"
#include "trm2F.h"
#include "logdecmp.h"

using namespace std;
using namespace GiNaC;

namespace xloops{

/*
	To decompose log, we need three things:
	1.	S function
	2.	z1sigma
	3.	z2sigma
*/

ex fn_S(int k, int l, int m, int n, const ex &sigma, const ex &z){
	ex return_S;

	return_S = mat_P[m][l][k]*sigma*pow(z,2) + (mat_E[m][l][k]+mat_Q[m][l][k]*sigma)*z - mat_msquare[k] + I*Rho;

	return return_S;
}

ex fn_z1sigma(int k, int l, int m, int n, const ex &sigma){
	ex return_z1sigma;

	ex a, b, c, delta;
	a = mat_P[m][l][k]*sigma;
	b = mat_E[m][l][k] + mat_Q[m][l][k]*sigma;
	c = -mat_msquare[k] + I*Rho;

	delta = b*b - 4.0*a*c;

	return_z1sigma = (-b - sqrt(delta)) / (2.0*a);

	return return_z1sigma;
}

ex fn_z2sigma(int k, int l, int m, int n, const ex &sigma){
	ex return_z2sigma;

	ex a, b, c, delta;
	a = mat_P[m][l][k]*sigma;
	b = mat_E[m][l][k] + mat_Q[m][l][k]*sigma;
	c = -mat_msquare[k] + I*Rho;

	delta = b*b - 4.0*a*c;

	return_z2sigma = (-b + sqrt(delta)) / (2.0*a);

	return return_z2sigma;
}

ex log_NPoint_82_plus(int k, int l, int m, int n, const ex &sigma, const ex &z){
	ex return_log_NPoint_82_plus;

	ex 	z1sigma = fn_z1sigma(k, l, m, n, sigma),
		z2sigma = fn_z2sigma(k, l, m, n, sigma),
		S = fn_S(k, l, m, n, sigma, z);

	return_log_NPoint_82_plus = log(S) + log(1.0 / (mat_P[m][l][k]*z + mat_Q[m][l][k]) );

	return_log_NPoint_82_plus += -2.0*I*Pi * 	                                         my_step(imag_part(mat_P[m][l][k]*z + mat_Q[m][l][k])) * my_step(-imag_part(S / (mat_P[m][l][k]*z + mat_Q[m][l][k])));

	return return_log_NPoint_82_plus;
}



ex log_NPoint_82_minus(int k, int l, int m, int n, const ex &sigma, const ex &z){
	ex return_log_NPoint_82_minus;

	ex 	z1sigma = fn_z1sigma(k, l, m, n, sigma),
		z2sigma = fn_z2sigma(k, l, m, n, sigma),
		S = fn_S(k, l, m, n, sigma, z);

	return_log_NPoint_82_minus = log(S) + log(-1.0 / (mat_P[m][l][k]*z + mat_Q[m][l][k]) );

	return_log_NPoint_82_minus += -2.0*I*Pi * my_step(imag_part( - (mat_P[m][l][k]*z + mat_Q[m][l][k]) )) * my_step(-imag_part(S / ( -(mat_P[m][l][k]*z + mat_Q[m][l][k]) )));

	return return_log_NPoint_82_minus;
}



ex log_Khiem_46_plus(int k, int l, int m, int n, const ex &sigma, const ex &z){
	ex log_plus = 0;

	ex 	z1sigma = fn_z1sigma(k, l, m, n, sigma),
		z2sigma = fn_z2sigma(k, l, m, n, sigma),
		S = fn_S(k, l, m, n, sigma, z);

	log_plus = log(mat_P[m][l][k]*sigma*z - mat_P[m][l][k]*sigma*z1sigma)
		+ log(z-z2sigma)
		- log(mat_P[m][l][k]*z + mat_Q[m][l][k]);

	log_plus += 2.0*Pi*I*(
	                       my_step(imag_part(mat_P[m][l][k]*sigma*z1sigma))
	                       * my_step(imag_part(z2sigma))
	                     )
		- 2.0*Pi*I*(
		             my_step(-imag_part(mat_Q[m][l][k]))
		             * my_step( - imag_part( S / (mat_P[m][l][k]*z + mat_Q[m][l][k]) ) )
		           );

	return log_plus;
}


ex log_Khiem_46_minus(int k, int l, int m, int n, const ex &sigma, const ex &z){
	ex log_minus = 0;

	ex 	z1sigma = fn_z1sigma(k, l, m, n, sigma),
		z2sigma = fn_z2sigma(k, l, m, n, sigma),
		S = fn_S(k, l, m, n, sigma, z);

	log_minus = log(-mat_P[m][l][k]*sigma*z + mat_P[m][l][k]*sigma*z1sigma)
		+ log(z-z2sigma)
		- log(mat_P[m][l][k]*z + mat_Q[m][l][k]);

	log_minus += -2.0*Pi*I*(
	                         my_step(imag_part(mat_P[m][l][k]*sigma*z1sigma))
	                         * my_step(-imag_part(z2sigma))
	                       )
		+ 2.0*Pi*I*(
		             my_step(imag_part(mat_Q[m][l][k]))
		             * my_step(-imag_part(
		                                   S / (mat_P[m][l][k]*z + mat_Q[m][l][k])
		                                 )
		                      )
		           );
	return log_minus;
}

ex log_original_plus(int k, int l, int m, int n, const ex &sigma, const ex &z){
	ex return_log_original_plus = 0;
	ex 	nom = fn_S(k, l, m, n, sigma, z),
		denom = mat_P[m][l][k]*z + mat_Q[m][l][k];
	if (denom.is_zero()){
		printf("\n\n Error in log_original, denominator is zero. (k=%d, l=%d, m=%d, n=%d)", k+1, l+1, m+1, n+1);
		cout << "\t\t z = " << z << endl << endl;
	}
	return_log_original_plus = log(nom / denom);

	return return_log_original_plus;
}

ex log_original_minus(int k, int l, int m, int n, const ex &sigma, const ex &z){
	ex return_log_original_minus = 0;
	ex 	nom = fn_S(k, l, m, n, sigma, z),
		denom = -( mat_P[m][l][k]*z + mat_Q[m][l][k] );
	if (denom.is_zero()){
		printf("\n\n Error in log_original, denominator is zero. (k=%d, l=%d, m=%d, n=%d)", k+1, l+1, m+1, n+1);
		cout << "\t\t z = " << z << endl << endl;
	}
	return_log_original_minus = log(nom / denom);

	return return_log_original_minus;
}
} // namespace xloops
