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
**						Khiem, Phan     (phanhongkhiem@gmail.com)
**						Quan, Phan      (anhquan.phanle@gmail.com)
********************************************************************************
**
**      Calculate NPoints terms
**
********************************************************************************
**
** Historial Log:
**      Date     Version Author   Description
**      ________________ _________________________________________
**      20090324 1.0     Quan PhanCreate this file
*******************************************************************************/

#include "ginac/ginac.h"
#include "extrm.h"
#include "my_fns.h"
#include "trmchk.h"
#include "lev_Npoint.h"

using namespace GiNaC;
namespace xloops{

ex fn_delta89(int m, int l, int k){
	ex delta89;
	// get the variables
	ex Qconj, E, P, mksquare;
	Qconj = mat_Q[m][l][k].conjugate();
	E = mat_E[m][l][k];
	P = mat_P[m][l][k];
	mksquare = mat_msquare[k];

	// calculate base on equ 89
	delta89 = 0;

	delta89 = pow(
	               imag_part(
	                          Qconj * E
	                          - P (mksquare - I*Rho)
	                        )
	               ,2 // square
	             );
	delta89 += 4.0 * imag_part( E*P ) * imag_part( Qconj*(mksquare - I*Rho) );

	return delta89;
}


ex fn_X189(int m, int l, int k){
	ex X189;
	// get the variables
	ex Qconj, E, P, mksquare;
	Qconj = mat_Q[m][l][k].conjugate();
	E = mat_E[m][l][k];
	P = mat_P[m][l][k];
	mksquare = mat_msquare[k];

	ex delta = fn_delta89(m, l, k);

	// calculate base on equ 89
	X189 = 0;

	X189 = -imag_part( Qconj*E - P*(mksquare - I*Rho) ) + sqrt(delta);
	X189 /= 2.0*imag_part(E*P);

	return X189;
}

ex fn_X289(int m, int l, int k){
	ex X289;
	// get the variables
	ex Qconj, E, P, mksquare;
	Qconj = mat_Q[m][l][k].conjugate();
	E = mat_E[m][l][k];
	P = mat_P[m][l][k];
	mksquare = mat_msquare[k];

	ex delta = fn_delta89(m, l, k);

	// calculate base on equ 89
	X289 = 0;

	X289 = -imag_part( Qconj*E - P*(mksquare - I*Rho) ) - sqrt(delta);
	X289 /= 2.0*imag_part(E*P);

	return X289;
}

ex fn_etaplus90 (int m, int l, int k, const ex &z){
	ex eta90 = 0;
	ex X1, X2;
	ex Qconj, E, P, mksquare;
	Qconj = mat_Q[m][l][k].conjugate();
	E = mat_E[m][l][k];
	P = mat_P[m][l][k];
	mksquare = mat_msquare[k];

	ex delta = fn_delta89(m, l, k);

	ex thetaMinusH;


	Qconj = mat_Q[m][l][k].conjugate();
	X1 = fn_X189(m, l, k);
	X2 = fn_X289(m, l, k);

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

	eta90 = -2.0 * I * Pi * my_step(imag_part(Qconj)) * thetaMinusH;
	return eta90;
}

ex fn_etaminus90 (int m, int l, int k, const ex &z){
	ex eta90 = 0;
	ex X1, X2;
	ex Qconj, E, P, mksquare;
	Qconj = mat_Q[m][l][k].conjugate();
	E = mat_E[m][l][k];
	P = mat_P[m][l][k];
	mksquare = mat_msquare[k];

	ex delta = fn_delta89(m, l, k);

	ex thetaPlusH;


	Qconj = mat_Q[m][l][k].conjugate();
	X1 = fn_X189(m, l, k);
	X2 = fn_X289(m, l, k);

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

	eta90 = -2.0 * I * Pi * my_step(imag_part(Qconj)) * thetaPlusH;
	return eta90;
}

}// Namespace xloops
