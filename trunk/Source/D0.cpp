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
**	Calculate D0, very important
**
********************************************************************************
**
** Historial Log:
**	Date		Version	Author		Description
**	_________	_______	_________	________________________________
**	20090214	1.0	Quan Phan	Create this file
*******************************************************************************/

#include "ginac/ginac.h"
#include "trm.h"
#include "my_fns.h"
#include "trmchk.h"
#include "lev1.h"
#include "lev2.h"
#include "lev3.h"
#include "lev4.h"
#include "lev5.h"
#include "ThetaG.h"
#include "RFunction.h"
#include "LogAG.h"

using namespace GiNaC;

namespace xloops{
	// need to calculate coff
ex fn_Coff(int n, int m, int l, int k){
	ex AC = fn_AC(l, k), Amlk = fn_A(m, l, k), Bmlk = fn_B(m, l, k), Anlk = fn_A(n, l, k), Bnlk = fn_B(n, l, k), beta = fn_beta(m, l, k), phi = fn_phi(m, l, k);
	ex coff = 0;

	coff = I*pow(Pi, 2) * (1.0/AC) * ( 1.0/(Bmlk*Anlk-Bnlk*Amlk) )
		* (1.0 - myfn_delta(AC)) * (1.0 - myfn_delta(Bmlk)) * abs(1.0 - beta*phi);

	return coff;
}

	// we need to calculate pos/ neg /ex term
ex fn_PosTerm(int n, int m, int l, int k){
	ex posTerm = 0;

			// Intermidate terms
	ex OPlus = fn_OPlus(n, m, l, k), T1 = fn_T1(n, m, l, k), T2 = fn_T2(n, m, l, k),
		f = fn_f(l, k), fminus = fn_fminus(l, k), g = fn_g(m, l, k), gminus = fn_gminus(m, l, k),
		beta = fn_beta(m, l, k), phi = fn_phi(m, l, k),
		F = fn_F(n, m, l, k), P = fn_P(m, l, k),
		z1beta = fn_z1beta(m, l, k), z2beta = fn_z2beta(m, l, k),
		z1phi = fn_z1phi(m, l, k), z2phi = fn_z2phi(m, l, k), Q = fn_Q(m, l, k);


			// Calculate it:
	posTerm = OPlus*RFunction(-T1, -T2) - f*g*LogAG((1.0-beta*phi)/beta, F/beta, -T1, -T2);
	posTerm += - f*gminus*LogAG(-(1.0-beta*phi)/beta, -F/beta, -T1, -T2) - (f*g+fminus*g) * LogAG(-P/beta, P*z1beta/beta, -T1, -T2);
	posTerm += -(f*g + fminus*g)*LogAG(1.0, -z2beta, -T1, -T2) + f*g*LogAG(-P*phi, P*phi*z1phi, -T1, -T2);
	posTerm += f*g*LogAG(1.00, -z2phi, -T1, -T2) * f*gminus*LogAG(P*phi, -P*phi*z1phi, -T1, -T2);
	posTerm += fminus*g*LogAG(1.0, -z2phi, -T1, -T2) + (fminus*g - f*gminus)*LogAG(P, Q, -T1, -T2);

	return posTerm;
}

ex fn_NegTerm(int n, int m, int l, int k){
	ex negTerm = 0;

			// Intermidate terms
	ex OMinus = fn_OMinus(n, m, l, k), T1 = fn_T1(n, m, l, k), T2 = fn_T2(n, m, l, k),
		f = fn_f(l, k), fminus = fn_fminus(l, k), g = fn_g(m, l, k), gminus = fn_gminus(m, l, k), beta = fn_beta(m, l, k), phi = fn_phi(m, l, k),
		F = fn_F(n, m, l, k), P = fn_P(m, l, k), z1beta = fn_z1beta(m, l, k), z2beta = fn_z2beta(m, l, k),
		z1phi = fn_z1phi(m, l, k), z2phi = fn_z2phi(m, l, k), Q = fn_Q(m, l, k);


			// Calculate it:
	negTerm = OMinus*RFunction(T1, T2) + (fminus*gminus + fminus*g)*LogAG(-(1.0-beta*phi)/beta, F/beta, T1, T2);
	negTerm += fminus*gminus*LogAG(P/beta, P*z1beta/beta, T1, T2) + fminus*gminus*LogAG(-1.0, -z2beta, T1, T2);
	negTerm += f*gminus*LogAG(-P/beta, -P*z1beta/beta, T1, T2) + f*gminus*LogAG(-1.0, -z2beta, T1, T2);
	negTerm += -(fminus*gminus + fminus*g)*LogAG(P*phi, P*phi*z1phi, T1, T2) - (fminus*gminus + fminus*g)*LogAG(-1.0, -z2phi, T1, T2);
	negTerm += (fminus*g - f*gminus) * LogAG(-P, Q, T1, T2);

	return negTerm;
}

ex fn_ExtraTerm(int n, int m, int l, int k){
	ex extraTerm = 0;
	ex f = fn_f(l, k), fminus = fn_fminus(l, k), g = fn_g(m, l, k), gminus = fn_gminus(m, l, k), ImQ = imag_part(fn_Q(m, l, k)), T1 = fn_T1(n, m, l, k), T2 = fn_T2(n, m, l, k), P = fn_P(m, l, k), E = fn_E(m, l, k), Q = fn_Q(m, l, k);
	ex A0 = imag_part(P*E);
	ex B0 = imag_part(E*Q.conjugate() - P*mat_msquare[k] + I*Rho*P);
	ex C0 = imag_part((-mat_msquare[k] + I*Rho)*Q.conjugate());
	extraTerm = 2.0*Pi*I*fminus*g*my_step(-ImQ)*ThetaG(A0, B0, C0, -T1, -T2);
	extraTerm += 2.0*Pi*I*f*gminus*my_step(ImQ)*ThetaG(-A0,-B0,-C0,-T1,-T2);

	return extraTerm;
}

ex D0(const ex& q_, const ex& m_, const ex& Rho_, const ex& Rho1_, const ex& Rho2_){
		// First thing first. init variable
	mat_msquare[0]=m_.op(0);	mat_msquare[1]=m_.op(1);	mat_msquare[2]=m_.op(2);	mat_msquare[3]=m_.op(3);
	mat_q[0][0]=q_.op(0);mat_q[0][1]=0;mat_q[0][2]=0;mat_q[0][3]=0;
	mat_q[1][0]=q_.op(1);mat_q[1][1]=q_.op(2);mat_q[1][2]=0;mat_q[1][3]=0;
	mat_q[2][0]=q_.op(3);mat_q[2][1]=q_.op(4);mat_q[2][2]=q_.op(5);mat_q[2][3]=0;
	Rho = Rho_; Rho1 = Rho1_; Rho2 = Rho2_;

		// Next to do: calculate the needed intermidate terms
	ex return_D0 = 0;

	for (int n=0; n<4; n++) for (int m=0; m<4; m++) for (int l=0; l<4; l++) for (int k=0; k<4; k++){
		if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k){
			ex coff = fn_C(m, l, k), posTerm = fn_PosTerm(n, m, l, k), negTerm = fn_NegTerm(n, m, l, k), extraTerm = fn_ExtraTerm(n, m, l, k);

			return_D0 += coff*(posTerm+negTerm+extraTerm);
		}
	}

	return return_D0;
}
}// Namespace xloops
