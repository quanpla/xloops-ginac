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

#include "trm.h"
#include "my_fns.h"
#include "trm2F.h"

#include "lev1.h"
#include "lev2.h"
#include "lev3.h"
#include "lev4.h"

#define ZRANGE 10
#define ZSTEP 1

using namespace GiNaC;
using namespace xloops;

/*
	To calculate equ 75, we need three things:
	1.	The constant
	2.	The positive integrand
	3.	The negative integrand
*/

ex fn_S(int k, int l, int m, int n, const ex &sigma, const ex &z){
	ex return_S;

	return_S = mat_P[m][l][k]*sigma*pow(z,2) + (mat_E[m][l][k]+mat_Q[m][l][k]*sigma)*z - mat_msquare[k] + I*Rho;

	return return_S;
}

ex fn_z1sigma(int k, int l, int m, int n, const ex &sigma){
	ex return_z1sigma;

	ex a, b, c, delta;
	a = 1;
	b = (mat_E[m][l][k] + mat_Q[m][l][k]*sigma) / (mat_P[m][l][k]*sigma);
	c = -mat_msquare[k] + I*Rho;

	delta = b*b - 4.0*a*c;

	return_z1sigma = (-b - sqrt(delta)) / (2.0*a);

	return return_z1sigma;
}

ex fn_z2sigma(int k, int l, int m, int n, const ex &sigma){
	ex return_z2sigma;

	ex a, b, c, delta;
	a = 1;
	b = (mat_E[m][l][k] + mat_Q[m][l][k]*sigma) / (mat_P[m][l][k]*sigma);
	c = -mat_msquare[k] + I*Rho;

	delta = b*b - 4.0*a*c;

	return_z2sigma = (-b + sqrt(delta)) / (2.0*a);

	return return_z2sigma;
}

ex log_plus_decompose(int k, int l, int m, int n, const ex &sigma, const ex &z){
	ex log_plus = 0;

	log_plus = log(mat_P[k][l][m]*sigma*z - mat_P[k][l][m]*sigma*fn_z1sigma(k, l, m, n, sigma))
		+ log(z-fn_z2sigma(k, l, m, n, sigma))
		- log(mat_P[k][l][m]*z + mat_Q[k][l][m]);

	log_plus += 2.0*Pi*I*(
	                       my_step(imag_part(mat_P[k][l][m]*sigma*fn_z1sigma(k, l, m, n, sigma)))
	                       * my_step(imag_part(fn_z2sigma(k, l, m, n, sigma)))
	                     )
		- 2.0*Pi*I*(
		             my_step(-imag_part(mat_Q[k][l][m]))
		             * my_step(imag_part(fn_S(k, l, m, n, sigma, z) / (mat_P[k][l][m]*z + mat_Q[k][l][m])))
		           );

	return log_plus;
}

ex log_minus_decompose(int k, int l, int m, int n, const ex &sigma, const ex &z){
	ex log_minus = 0;

	log_minus = log(-mat_P[k][l][m]*sigma*z + mat_P[k][l][m]*sigma*fn_z1sigma(k, l, m, n, sigma))
		+ log(z-fn_z2sigma(k, l, m, n, sigma))
		- log(mat_P[k][l][m]*z + mat_Q[k][l][m]);

	log_minus += -2.0*Pi*I*(
	                       my_step(imag_part(mat_P[k][l][m]*sigma*fn_z1sigma(k, l, m, n, sigma)))
	                       * my_step(-imag_part(fn_z2sigma(k, l, m, n, sigma)))
	                     )
		+ 2.0*Pi*I*(
		             my_step(imag_part(mat_Q[k][l][m]))
		             * my_step(-imag_part(
		                                     fn_S(k, l, m, n, sigma, z)
		                                     /(mat_P[k][l][m]*z + mat_Q[k][l][m])
		                                   )
		                        )
		           );
	return log_minus;
}

int main(int argc, char * argv[]){
	double p10, p20, p21, p30, p31, p32, m1s, m2s, m3s, m4s;

	p10 = atof(argv[1]);
	p20 = atof(argv[2]);
	p21 = atof(argv[3]);
	p30 = atof(argv[4]);
	p31 = atof(argv[5]);
	p32 = atof(argv[6]);
	m1s = atof(argv[7]);
	m2s = atof(argv[8]);
	m3s = atof(argv[9]);
	m4s = atof(argv[10]);

	ex p, q, msquare;

	msquare = lst(m1s, m2s, m3s, m4s);
	p = lst(p10, p20, p21, p30, p31, p32);
        // convert p --> q
	q = p2q(p);

	Rho = 1e-20; Rho1 = 0; Rho2 = 0;

	// from here is xloops-ginac calculation
	mat_msquare[0]=m1s; mat_msquare[1]=m2s; mat_msquare[2]=m3s; mat_msquare[3]=m4s;
	mat_q[0][0]=q.op(0);mat_q[0][1]=0;mat_q[0][2]=0;mat_q[0][3]=0;
	mat_q[1][0]=q.op(1);mat_q[1][1]=q.op(2);mat_q[1][2]=0;mat_q[1][3]=0;
	mat_q[2][0]=q.op(3);mat_q[2][1]=q.op(4);mat_q[2][2]=q.op(5);mat_q[2][3]=0;
	input_Rho = Rho; input_Rho1 = Rho1; input_Rho2 = Rho2;

	xloopsGiNaC_calc_lev1();
	xloopsGiNaC_calc_lev2();
	xloopsGiNaC_calc_lev3();
	xloopsGiNaC_calc_lev4();

	int k, l, m, n;
	ex z = symbol("z");
	n = 0; // don't care n for now
	for(k=0; k<4; k++) for(l=0; l<4; l++) for(m = 0; m<4; m++){
		if(m!=l && m!=k && l!=k)
		{
			printf( "(m,l,k) = (%d,%d,%d) \n\t\t ", m,l,k);

			cout << mat_beta[m][l][k] << "\t" << mat_phi[m][l][k] << "\n" ;

			ex logbeta_plus_1, logphi_plus_1, logbeta_plus_2, logphi_plus_2; // log_plus_1 = the compact form of log, log_plus_2 = the decompose form of log
			ex logbeta_minus_1, logphi_minus_1, logbeta_minus_2, logphi_minus_2;

			logbeta_plus_1 = log(fn_S(k, l, m, n, mat_beta[m][l][k], z)
			                     /(mat_P[m][l][k]*z + mat_Q[m][l][k])
			                    );
			logbeta_plus_2 = log_plus_decompose(k, l, m, n, mat_beta[m][l][k], z);
			logphi_plus_1 = log(fn_S(k, l, m, n, mat_phi[m][l][k], z)
			                     /(mat_P[m][l][k]*z + mat_Q[m][l][k])
			                    );
			logphi_plus_2 = log_plus_decompose(k, l, m, n, mat_phi[m][l][k], z);

			logbeta_minus_1 = log(-fn_S(k, l, m, n, mat_beta[m][l][k], z)
			                     /(mat_P[m][l][k]*z + mat_Q[m][l][k])
			                    );
			logbeta_minus_2 = log_minus_decompose(k, l, m, n, mat_beta[m][l][k], z);
			logphi_minus_1 = log(-fn_S(k, l, m, n, mat_phi[m][l][k], z)
			                    /(mat_P[m][l][k]*z + mat_Q[m][l][k])
			                   );
			logphi_minus_2 = log_minus_decompose(k, l, m, n, mat_phi[m][l][k], z);

			for(double dz=-ZRANGE; dz<=ZRANGE; dz+=ZSTEP){
				// z, logbeta_1+z, logbeta_2+z, dif, ...
				printf("\n%e", dz);
				ex val = logbeta_plus_1.subs(z==dz);
				cout << "\t" << val.evalf();
				val = logbeta_plus_2.subs(z==dz);
				cout << "\t" << val.evalf();
				val = (logbeta_plus_1-logbeta_plus_2).subs(z==dz);
				cout << "\t" << val.evalf();

				val = logbeta_minus_1.subs(z==dz);
				cout << "\t" << val.evalf();
				val = logbeta_minus_2.subs(z==dz);
				cout << "\t" << val.evalf();
				val = (logbeta_minus_1-logbeta_plus_2).subs(z==dz);
				cout << "\t" << val.evalf();

				val = logphi_plus_1.subs(z==dz);
				cout << "\t" << val.evalf();
				val = logphi_plus_2.subs(z==dz);
				cout << "\t" << val.evalf();
				val = (logphi_plus_1-logphi_plus_2).subs(z==dz);
				cout << "\t" << val.evalf();

				val = logphi_minus_1.subs(z==dz);
				cout << "\t" << val.evalf();
				val = logphi_minus_2.subs(z==dz);
				cout << "\t" << val.evalf();
				val = (logphi_minus_1-logphi_plus_2).subs(z==dz);
				cout << "\t" << val.evalf();
			}

		}
	}


	return EXIT_SUCCESS;
}
