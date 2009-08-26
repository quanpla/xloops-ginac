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
ex Jacobian(int n, int m, int l, int k){
	ex beta_mlk = fn_beta(m, l, k), phi_mlk = fn_phi(m, l, k), AC_lk = fn_AC(l, k),
		A_mlk = fn_A(m, l, k), A_nlk = fn_A(n, l, k),
		B_mlk = fn_B(m, l, k), B_nlk = fn_B(n, l, k),
		P_mlk = fn_P(m, l, k);
	
	return abs(1.0 - beta_mlk*phi_mlk) 
		* (1.0 - myfn_delta(AC_lk)) 
		* (1.0 - myfn_delta(B_mlk)) 
		* (1.0/AC_lk)
		* (1.0/( B_mlk*A_nlk-B_nlk*A_mlk ))
		* (-1.0/P_mlk);
}

ex D0_part(int n, int m, int l, int k){
	/* A part of D0 with index n, m, l, k*/
	ex term;
	
	ex 	OPlus = fn_OPlus(n, m, l, k), OMinus = fn_OMinus(n, m, l, k),
		T1 = fn_T1(n, m, l, k),	T2 = fn_T2(n, m, l, k),
		f = fn_f(l, k), fminus = fn_fminus(l, k), g = fn_g(m, l, k), gminus = fn_gminus(m, l, k),
		beta = fn_beta(m, l, k), phi = fn_phi(m, l, k),
		z1beta = fn_z1beta(m, l, k), z2beta = fn_z2beta(m, l, k), z1phi = fn_z1phi(m, l, k), z2phi = fn_z2phi(m, l, k),
		Q_im = fn_Q_im(m, l, k),
		A0 = fn_A0(m, l, k), B0 = fn_B0(m, l, k), C0 = fn_C0(m, l, k),
		P = fn_P(m, l, k), Q = fn_Q(m, l, k),
		F = fn_F(n, m, l, k);
	
	//	1.	No "Log" scalar
	term = 	OPlus*RFunction(-T1, -T2)
		+ 	OMinus*RFunction(T1, T2);
	
	// 	2.	Scalars with F_nmlk
	term += -	f*g*LogAG( /*1st*/ ( 1.0 - beta*phi )/beta, /*2nd*/ F/beta, /*3rd*/ -T1, /*4th*/ -T2)
		+	(fminus*gminus + fminus*g) * LogAG(/*1st*/ -( 1.0 - beta*phi )/beta, /*2nd*/ F/beta, /*3rd*/ T1, /*4th*/ T2)
		-	f*gminus * LogAG(/*1st*/ -( 1.0 - beta*phi )/beta, /*2nd*/ -F/beta, /*3rd*/ -T1, /*4th*/ -T2);
	
	//	3.	Scalars with log but no F_nmlk
	term += -(f*g + fminus*g) * LogAG( -P/beta, P*z1beta/beta, -T1, -T2)
		-	(f*g + fminus*g) * LogAG(1.0, -z2beta, -T1, -T2)
		+	f*g * LogAG(-P*phi, P*phi*z1phi, -T1, -T2)
		+	f*g * LogAG(1.0, -z2phi, -T1, -T2)
		+	f*gminus * LogAG(P*phi, -P*phi*z1phi, -T1, -T2)
		+	f*gminus * LogAG(1.0, -z2phi, -T1, -T2)
		+	(fminus*g - f*gminus) * LogAG(P, Q, -T1, -T2)
		+	fminus*gminus * LogAG(P/beta, P*z1beta/beta, T1, T2)
		+	fminus*gminus * LogAG(-1.0, -z2beta, T1, T2)
		+	f*gminus * LogAG(-P/beta, -P*z1beta/beta, T1, T2)
		+ 	f*gminus * LogAG(-1.0, -z2beta, T1, T2)
		-	(fminus*gminus + fminus*g) * LogAG(P*phi, P*phi*z1phi, T1, T2)
		-	(fminus*gminus + fminus*g) * LogAG(-1.0, -z2phi, T1, T2)
		+	(fminus*g - f*gminus) * LogAG(-P, Q, T1, T2);
	
	//	4.	For M complex
	term += my_step(Q_im) * fminus*g + my_step(Q_im)*f*g* ThetaGC(-A0, -B0, -C0, -T1, -T2) * 2.0 * Pi * I;
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
			ex jacobi = Jacobian(n, m, l, k);
			return_D0 += jacobi*D0_part(n, m, l, k);
		}
	}

	return return_D0;
}
}// Namespace xloops
