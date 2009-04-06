/*******************************************************************************
**
** xloop-ginacs Project
** 1Loop4Pt implementation
**
**
** HCMUNS,
** June 2008
**
**
** Author(s): Son, D. H. (sondo01@gmail.com)
**						Khiem, Phan (phanhongkhiem@gmail.com)
**						Quan, Phan (anhquan.phanle@gmail.com)
********************************************************************************
**
** Calculate NPoints terms
**
********************************************************************************
**
** Historial Log:
** Date Version Author Description
** ________________ _________________________________________
** 20090324 1.0 Quan PhanCreate this file
*******************************************************************************/

#include "ginac/ginac.h"
#include "extrm.h"
#include "my_fns.h"
#include "trmchk.h"
#include "lev_Npoint.h"

using namespace GiNaC;
namespace xloops{

ex fn_T1(int n, int m, int n, int k){
	// T1, T2, ... for GZ function
	ex T1 = 0;
	const ex Q = mat_Q[m][l][k], ex F = mat_F[n][m][l][k], ex P = mat_P[m][l][k], ex beta = mat_beta[m][l][k],
		ex E = mat_E[m][l][k], m2k = mat_msquare[k];

	T1 = -(Q + F*P - beta*E) + sqrt( pow(Q + F*P - beta*E, 2) - 4.0*P*(Q*F + beta*(m2k - I*Rho) ) );
	T1/= 2.0*P;
	return T1;
}

ex fn_T2(int n, int m, int n, int k){
	// T1, T2, ... for GZ function
	ex T2 = 0;
	const ex Q = mat_Q[m][l][k], ex F = mat_F[n][m][l][k], ex P = mat_P[m][l][k], ex beta = mat_beta[m][l][k],
		ex E = mat_E[m][l][k], m2k = mat_msquare[k];

	T2 = -(Q + F*P - beta*E) - sqrt( pow(Q + F*P - beta*E, 2) - 4.*P*(Q*F + beta*(m2k - I*Rho) ) );
	T2/= 2.*P;
	return T2;
}

ex fn_T3(int n, int m, int n, int k){
	// T1, T2, ... for GZ function
	ex T3 = 0;
	const ex Q = mat_Q[m][l][k], ex F = mat_F[n][m][l][k], ex P = mat_P[m][l][k], ex beta = mat_beta[m][l][k],
		ex E = mat_E[m][l][k], m2k = mat_msquare[k];

	T3 = (Q + F*P - beta*E) + sqrt( pow(Q + F*P - beta*E, 2) - 4.0*P*(Q*F + beta*(m2k - I*Rho) ) );
	T3/= 2.0*P;
	return T3;
}

ex fn_T4(int n, int m, int n, int k){
	// T1, T2, ... for GZ function
	ex T4 = 0;
	const ex Q = mat_Q[m][l][k], ex F = mat_F[n][m][l][k], ex P = mat_P[m][l][k], ex beta = mat_beta[m][l][k],
		ex E = mat_E[m][l][k], m2k = mat_msquare[k];

	T4 = (Q + F*P - beta*E) - sqrt( pow(Q + F*P - beta*E, 2) - 4.0*P*(Q*F + beta*(m2k - I*Rho) ) );
	T4/= 2.0*P;
	return T4;
}

ex fn_R_m1e(const ex &epsilon, const ex &x, const ex &y, const ex &z){
 /*
  Function R_(-1-e) [1, 1, eps; x, y, z]
 */
	ex R_m1e = 0;
	ex logx = log(x), logy = log(y);

 /* line 1 */
	R_m1e = (logx-logy)/(x-y) + epsilon*(logx-logy)/(x-y);

 /* orther lines */
	R_m1e += ( epsilon/(x-y) )/
		(
		  pow(logy, 2)/2.0 - pow(logx, 2)/2.0 + Li2(1.0 - z/y) - Li2(1.0 - z/x)
		  + logy*( fn_eta(z-y, 1.0/(1.0-y) ) - fn_eta(z-y, -1.0/y) )
		  - logx*( fn_eta(z-x, 1.0/(1.0-x) ) - fn_eta(z-x, -1.0/x) )
		  + log(1.0-z/y)*fn_eta(z, 1.0/y) - log(1.0-z/x)*fn_eta(z, 1.0/x)
		);

	return R_m1e;
}

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

// Z1 and Z2 in calculating log(S/ (Pz+Q) )
ex fn_Z1(int n, int m, int l, int k, const ex &sigma){
	ex Z1 = 0;
	ex 	E = mat_E[m][l][k],
		P = mat_P[m][l][k],
		Q = mat_Q[m][l][k],
		m2k = mat_msquare[k];

	Z1 = -(E + Q*sigma)	+ sqrt(
	                   	        pow(E + Q*sigma, 2) + 4.0*P*sigma*(m2k - I*Rho)
	                   	      );
	Z1 /= 2.0*P*sigma;

	return Z1;
}
ex fn_Z2(int n, int m, int l, int k, const ex &sigma){
	ex Z2 = 0;
	ex 	E = mat_E[m][l][k],
		P = mat_P[m][l][k],
		Q = mat_Q[m][l][k],
		m2k = mat_msquare[k];

	Z2 = -(E + Q*sigma)	- sqrt(
	                   	        pow(E + Q*sigma, 2) + 4.0*P*sigma*(m2k - I*Rho)
	                   	      );
	Z2 /= 2.0*P*sigma;

	return Z2;
}

ex fn_delta90(int m, int l, int k){
	ex delta90;
	// get the variables
	ex Qconj, E, P, mksquare;
	Qconj = mat_Q[m][l][k].conjugate();
	E = mat_E[m][l][k];
	P = mat_P[m][l][k];
	mksquare = mat_msquare[k];

	// calculate base on equ 90
	delta90 = 0;

	delta90 = pow(
	               imag_part(
	                          Qconj * E
	                          - P (mksquare - I*Rho)
	                        )
	               ,2 // square
	             );
	delta90 += 4.0 * imag_part( E*P ) * imag_part( Qconj*(mksquare - I*Rho) );

	return delta90;
}


ex fn_X190(int m, int l, int k){
	ex X190;
	// get the variables
	ex Qconj, E, P, mksquare;
	Qconj = mat_Q[m][l][k].conjugate();
	E = mat_E[m][l][k];
	P = mat_P[m][l][k];
	mksquare = mat_msquare[k];

	ex delta = fn_delta90(m, l, k);

	// calculate base on equ 90
	X190 = 0;

	X190 = -imag_part( Qconj*E - P*(mksquare - I*Rho) ) + sqrt(delta);
	X190 /= 2.0*imag_part(E*P);

	return X190;
}

ex fn_X290(int m, int l, int k){
	ex X290;
	// get the variables
	ex Qconj, E, P, mksquare;
	Qconj = mat_Q[m][l][k].conjugate();
	E = mat_E[m][l][k];
	P = mat_P[m][l][k];
	mksquare = mat_msquare[k];

	ex delta = fn_delta90(m, l, k);

	// calculate base on equ 90
	X290 = 0;

	X290 = -imag_part( Qconj*E - P*(mksquare - I*Rho) ) - sqrt(delta);
	X290 /= 2.0*imag_part(E*P);

	return X290;
}

ex fn_etaplus100 (int m, int l, int k, const ex &z){
	ex eta100 = 0;
	ex X1, X2;
	ex Qconj, E, P, mksquare;
	Qconj = mat_Q[m][l][k].conjugate();
	E = mat_E[m][l][k];
	P = mat_P[m][l][k];
	mksquare = mat_msquare[k];

	ex delta = fn_delta90(m, l, k);

	ex thetaMinusH;


	Qconj = mat_Q[m][l][k].conjugate();
	X1 = fn_X190(m, l, k);
	X2 = fn_X290(m, l, k);

	// now calculate thetaMinusH (91) to (97)
	thetaMinusH = 0;
	if( imag_part(E*P).info(info_flags::positive) ){
		if ( X1 <= z && z <= X2 && delta >= 0.0)
			thetaMinusH = 1; // equ 91
	}
	else{
		if(imag_part(E*P).info(info_flags::negative)){
			if ( z <= X1 && X2 <= z && delta.info(info_flags::positive))
				thetaMinusH = 1; // equ 93
			else{
				if (delta.info(info_flags::negative)){
					thetaMinusH = 1; // equ 93
				}
			}
		}
		else{
			if(imag_part(Qconj*E - P*(mksquare - I*Rho) ) != 0){
				ex compareCondition;
				compareCondition = imag_part(Qconj* (mksquare - I*Rho) ) / (Qconj*E - P*(mksquare - I*Rho) )
					if (z <= compareCondition)
						thetaMinusH = 1; // equ 95
				else
					thetaMinusH my_step(imag_part(Qconj*(mksquare - I*Rho)) ); // equ 97
			}
		}
	}

	eta100 = -2.0 * I * Pi * my_step(imag_part(Qconj)) * thetaMinusH;
	return eta100;
}

ex fn_etaminus100 (int m, int l, int k, const ex &z){
	ex eta100 = 0;
	ex X1, X2;
	ex Qconj, E, P, mksquare;
	Qconj = mat_Q[m][l][k].conjugate();
	E = mat_E[m][l][k];
	P = mat_P[m][l][k];
	mksquare = mat_msquare[k];

	ex delta = fn_delta90(m, l, k);

	ex thetaPlusH;


	Qconj = mat_Q[m][l][k].conjugate();
	X1 = fn_X190(m, l, k);
	X2 = fn_X290(m, l, k);

	// now calculate thetaPlusH (91) to (97)
	thetaPlusH = 0;
	if( imag_part(E*P).info(info_flags::positive) ){
		if ( X1 <= z && z <= X2 && delta >= 0.0)
			thetaPlusH = 1; // equ 91
	}
	else{
		if(imag_part(E*P).info(info_flags::negative)){
			if ( z <= X1 && X2 <= z && delta.info(info_flags::positive))
				thetaPlusH = 1; // equ 93
			else{
				if (delta.info(info_flags::negative)){
					thetaPlusH = 1; // equ 93
				}
			}
		}
		else{
			if(imag_part(Qconj*E - P*(mksquare - I*Rho) ) != 0){
				ex compareCondition;
				compareCondition = imag_part(Qconj* (mksquare - I*Rho) ) / (Qconj*E - P*(mksquare - I*Rho) )
					if (z <= compareCondition)
						thetaPlusH = 1; // equ 95
				else
					thetaPlusH my_step(imag_part(Qconj*(mksquare - I*Rho)) ); // equ 97
			}
		}
	}

	eta100 = -2.0 * I * Pi * my_step(imag_part(Qconj)) * thetaPlusH;
	return eta100;
}

ex fn_roundplus(int n, int m, int l, int k){
	ex roundplus = 0;

	for(int ki = 0; ki<4; ki++)for(int li = 0; li<4; li++)for(int mi = 0; mi<4; mi++){
		if(ki!=li && ki!=mi && li!=mi){
			roundplus += /*nom*/ (1.0 - myfn_delta(mat_AC[l][k]) )  * (1.0 - myfn_delta(mat_B[m][l][k]))
				/*denom*/ / mat_AC[l][k]*( mat_B[m][l][k]*mat_A[n][l][k] - mat_B[n][l][k]*mat_A[m][l][k]);
		}
	}

	roundplus = I*Pow(Pi, 2);

	return roundplus;
}

ex fn_Eta(int m, int l, int k, const ex &sigma){
	ex Eta = 0;

	const ex P = mat_P[k][l][m],
		Z1 = fn_Z1(n, m, l, k, sigma),
		Z2 = fn_Z2(n, m, l, k, sigma);

	Eta = 2.0*I*Pi*my_step( imag_part(P*sigma*Z1) ) * my_step( imag_part(Z2) );

	return Eta;
}

ex fn_B(const ex &x, const ex &y){
	// B(1+x, y) = 1+x for now
	return x;
}

ex fn_varLplus(const ex &a, const ex &b, const ex & epsilon){
	ex varLplus = 0;

	const ex T1 = fn_T1(n, m, l, k), T2 = fn_T2(n, m, l, k), T3 = fn_T3(n, m, l, k), T4 = fn_T4(n, m, l, k);
	//fn_R_m1e(const ex &epsilon, const ex &x, const ex &y, const ex &z){
	ex eps = possymbol("epsilon");

	varLplus = 1.0/eps;
	varLplus *= ( fn_GZ(T3, T4) + fn_B(1.0 + eps, 1.0)*pow(a,-eps)*fn_R_m1e(eps, -T1, -T2, b/a) );

	varLplus = varLplus.subs(eps==epsilon);
	return varLplus;
}

ex fn_varLminus(const ex &a, const ex &b, const ex & epsilon){
	ex varLminus = 0;

	const ex T1 = fn_T1(n, m, l, k), T2 = fn_T2(n, m, l, k), T3 = fn_T3(n, m, l, k), T4 = fn_T4(n, m, l, k);
	//fn_R_m1e(const ex &epsilon, const ex &x, const ex &y, const ex &z){
	ex eps = possymbol("epsilon");

	varLminus = 1.0/eps;
	varLminus *= ( fn_GZ(T1, T2) + fn_B(1.0 + eps, 1.0)*pow(-a,-eps)*fn_R_m1e(eps, -T3, -T4, -b/a) );

	varLminus = varLminus.subs(eps==epsilon);
	return varLminus;
}

}// Namespace xloops
