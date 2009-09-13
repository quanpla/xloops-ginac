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
void trm2F();
void fn2F();
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
	//printf("1. No log scalar.\n");
	//	1.	No "Log" scalar
	term = 	OPlus*RFunction(-T1, -T2)
		+ 	OMinus*RFunction(T1, T2);
	
	//printf("2. Scalars with F_nmlk.\n");
	// 	2.	Scalars with F_nmlk
	term += -	f*g*LogAG( /*1st*/ ( 1.0 - beta*phi )/beta, /*2nd*/ F/beta, /*3rd*/ -T1, /*4th*/ -T2)
		+	(fminus*gminus + fminus*g) * LogAG(/*1st*/ -( 1.0 - beta*phi )/beta, /*2nd*/ F/beta, /*3rd*/ T1, /*4th*/ T2)
		-	f*gminus * LogAG(/*1st*/ -( 1.0 - beta*phi )/beta, /*2nd*/ -F/beta, /*3rd*/ -T1, /*4th*/ -T2);
	
	//printf("3. Scalars with log but no F_nmlk.\n");
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
	
	//printf("4. For M complex.\n");
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
	//cout << endl << "Start to calculate" << endl;
		// Next to do: calculate the needed intermidate terms
	ex return_D0 = 0;

	// init intermidiate terms
	lev4Calc();
	//for (int n=0; n<4; n++) for (int m=0; m<4; m++) for (int l=0; l<4; l++) for (int k=0; k<4; k++){
	//	if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k){
	//		printf("1.    %d\t%d\t%d\t%d\n", n, m, l, k);
	//		ex jacobi = Jacobian(n, m, l, k);
	//		printf("2.    %d\t%d\t%d\t%d\n", n, m, l, k);
	//		return_D0 += jacobi*D0_part(n, m, l, k);
	//	}
	//}
	//cout << "*****************************************************\nExport for view";
	fn2F();
	//cout << endl << endl;
	//cout << "End of Export\n";
	return return_D0;
}
void trm2F(){
	int n, m, l, k;
		cout <<"Quan's Code" << endl << "q_ij" << endl;
		cout << evalf(mat_q[0][0]) << "\t" << evalf(mat_q[0][1]) << "\t" << evalf(mat_q[0][2]) << "\t" << evalf(mat_q[0][3]) << "\t" << endl;
		cout << evalf(mat_q[1][0]) << "\t" << evalf(mat_q[1][1]) << "\t" << evalf(mat_q[1][2]) << "\t" << evalf(mat_q[1][3]) << "\t" << endl;
		cout << evalf(mat_q[2][0]) << "\t" << evalf(mat_q[2][1]) << "\t" << evalf(mat_q[2][2]) << "\t" << evalf(mat_q[2][3]) << "\t" << endl;
		cout << evalf(mat_q[3][0]) << "\t" << evalf(mat_q[3][1]) << "\t" << evalf(mat_q[3][2]) << "\t" << evalf(mat_q[3][3]) << "\t" << endl;
		cout <<"m_i square" << endl;
		cout << evalf(mat_msquare[0]) << "\t" << evalf(mat_msquare[1]) << "\t" << evalf(mat_msquare[2]) << "\t" << evalf(mat_msquare[3]) << "\t" << endl;

//		for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
//			cout << "a_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_a[l][k])) << "\t" << imag_part(evalf(mat_a[l][k])) << endl;
//		for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
//			cout << "b_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_b[l][k])) << "\t" << imag_part(evalf(mat_b[l][k])) << endl;
//		for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
//			cout << "c_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_c[l][k])) << "\t" << imag_part(evalf(mat_c[l][k])) << endl;
//		for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
//			cout << "d_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_d[l][k])) << "\t" << imag_part(evalf(mat_d[l][k])) << endl;
//		for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
//			cout << "AC_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_AC[l][k])) << "\t" << imag_part(evalf(mat_AC[l][k])) << endl;
//		for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
//			cout << "alpha_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_alpha[l][k])) << "\t" << imag_part(evalf(mat_alpha[l][k])) << endl;
//		for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
//			cout << "f_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_f[l][k])) << "\t" << imag_part(evalf(mat_f[l][k])) << endl;
//		for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
//			cout << "fminus_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_fminus[l][k])) << "\t" << imag_part(evalf(mat_fminus[l][k])) << endl;
//		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
//			cout << "A_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_A[m][l][k])) << "\t" << imag_part(evalf(mat_A[m][l][k])) << endl;
//		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
//			cout << "B_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_B[m][l][k])) << "\t" << imag_part(evalf(mat_B[m][l][k])) << endl;
//		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
//			cout << "C_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_C[m][l][k])) << "\t" << imag_part(evalf(mat_C[m][l][k])) << endl;
//		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
//			cout << "D_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_D[m][l][k])) << "\t" << imag_part(evalf(mat_D[m][l][k])) << endl;
//		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
//			cout << "beta_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_beta[m][l][k])) << "\t" << imag_part(evalf(mat_beta[m][l][k])) << endl;
//		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
//			cout << "phi_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_phi[m][l][k])) << "\t" << imag_part(evalf(mat_phi[m][l][k])) << endl;
//		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
//			cout << "Q_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_Q[m][l][k])) << "\t" << imag_part(evalf(mat_Q[m][l][k])) << endl;
//		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
//			cout << "P_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_P[m][l][k])) << "\t" << imag_part(evalf(mat_P[m][l][k])) << endl;
//		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
//			cout << "E_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_E[m][l][k])) << "\t" << imag_part(evalf(mat_E[m][l][k])) << endl;
//		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
//			cout << "g_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_g[m][l][k])) << "\t" << imag_part(evalf(mat_g[m][l][k])) << endl;	
//		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
//			cout << "gminus_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_gminus[m][l][k])) << "\t" << imag_part(evalf(mat_gminus[m][l][k])) << endl;
//		for(n=0; n<4; n++) for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k)
//			cout << "F_" << "\t" << n+1 << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_F[n][m][l][k])) << "\t" << imag_part(evalf(mat_F[n][m][l][k])) << endl;
//		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
//			cout << "z1beta_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_z1beta[m][l][k])) << "\t" << imag_part(evalf(mat_z1beta[m][l][k])) << endl;
//		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
//			cout << "z2beta_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_z2beta[m][l][k])) << "\t" << imag_part(evalf(mat_z2beta[m][l][k])) << endl;
//		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
//			cout << "z1phi_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_z1phi[m][l][k])) << "\t" << imag_part(evalf(mat_z1phi[m][l][k])) << endl;
//		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
//			cout << "z2phi_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_z2phi[m][l][k])) << "\t" << imag_part(evalf(mat_z2phi[m][l][k])) << endl;
//		for(n=0; n<4; n++) for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k)
//			cout << "T1_" << "\t" << n+1 << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_T1[n][m][l][k])) << "\t" << imag_part(evalf(mat_T1[n][m][l][k])) << endl;
//		for(n=0; n<4; n++) for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k)
//			cout << "T2_" << "\t" << n+1 << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_T2[n][m][l][k])) << "\t" << imag_part(evalf(mat_T2[n][m][l][k])) << endl;
//		for(n=0; n<4; n++) for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k)
//			cout << "OMinus_" << "\t" << n+1 << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_OMinus[n][m][l][k])) << "\t" << imag_part(evalf(mat_OMinus[n][m][l][k])) << endl;
//		for(n=0; n<4; n++) for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k)
//			cout << "OPlus_" << "\t" << n+1 << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_OPlus[n][m][l][k])) << "\t" << imag_part(evalf(mat_OPlus[n][m][l][k])) << endl;
//		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
//			cout << "A0_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_A0[m][l][k])) << "\t" << imag_part(evalf(mat_A0[m][l][k])) << endl;
//		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
//			cout << "B0_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_B0[m][l][k])) << "\t" << imag_part(evalf(mat_B0[m][l][k])) << endl;
//		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
//			cout << "C0_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(mat_C0[m][l][k])) << "\t" << imag_part(evalf(mat_C0[m][l][k])) << endl;
//
		for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
			cout << "a_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_a[l][k])) << "\t" << imag_part((mat_a[l][k])) << endl;
		for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
			cout << "b_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_b[l][k])) << "\t" << imag_part((mat_b[l][k])) << endl;
		for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
			cout << "c_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_c[l][k])) << "\t" << imag_part((mat_c[l][k])) << endl;
		for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
			cout << "d_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_d[l][k])) << "\t" << imag_part((mat_d[l][k])) << endl;
		for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
			cout << "AC_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_AC[l][k])) << "\t" << imag_part((mat_AC[l][k])) << endl;
		for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
			cout << "alpha_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_alpha[l][k])) << "\t" << imag_part((mat_alpha[l][k])) << endl;
		for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
			cout << "f_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_f[l][k])) << "\t" << imag_part((mat_f[l][k])) << endl;
		for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
			cout << "fminus_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_fminus[l][k])) << "\t" << imag_part((mat_fminus[l][k])) << endl;
		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
			cout << "A_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_A[m][l][k])) << "\t" << imag_part((mat_A[m][l][k])) << endl;
		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
			cout << "B_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_B[m][l][k])) << "\t" << imag_part((mat_B[m][l][k])) << endl;
		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
			cout << "C_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_C[m][l][k])) << "\t" << imag_part((mat_C[m][l][k])) << endl;
		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
			cout << "D_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_D[m][l][k])) << "\t" << imag_part((mat_D[m][l][k])) << endl;
/*		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
			cout << "beta_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_beta[m][l][k])) << "\t" << imag_part((mat_beta[m][l][k])) << endl;
		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
			cout << "phi_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_phi[m][l][k])) << "\t" << imag_part((mat_phi[m][l][k])) << endl;
		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
			cout << "Q_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_Q[m][l][k])) << "\t" << imag_part((mat_Q[m][l][k])) << endl;
		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
			cout << "P_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_P[m][l][k])) << "\t" << imag_part((mat_P[m][l][k])) << endl;
		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
			cout << "E_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_E[m][l][k])) << "\t" << imag_part((mat_E[m][l][k])) << endl;
		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
			cout << "g_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_g[m][l][k])) << "\t" << imag_part((mat_g[m][l][k])) << endl;	
		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
			cout << "gminus_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_gminus[m][l][k])) << "\t" << imag_part((mat_gminus[m][l][k])) << endl;
		for(n=0; n<4; n++) for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k)
			cout << "F_" << "\t" << n+1 << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_F[n][m][l][k])) << "\t" << imag_part((mat_F[n][m][l][k])) << endl;
		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
			cout << "z1beta_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_z1beta[m][l][k])) << "\t" << imag_part((mat_z1beta[m][l][k])) << endl;
		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
			cout << "z2beta_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_z2beta[m][l][k])) << "\t" << imag_part((mat_z2beta[m][l][k])) << endl;
		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
			cout << "z1phi_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_z1phi[m][l][k])) << "\t" << imag_part((mat_z1phi[m][l][k])) << endl;
		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
			cout << "z2phi_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_z2phi[m][l][k])) << "\t" << imag_part((mat_z2phi[m][l][k])) << endl;
		for(n=0; n<4; n++) for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k)
			cout << "T1_" << "\t" << n+1 << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_T1[n][m][l][k])) << "\t" << imag_part((mat_T1[n][m][l][k])) << endl;
		for(n=0; n<4; n++) for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k)
			cout << "T2_" << "\t" << n+1 << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_T2[n][m][l][k])) << "\t" << imag_part((mat_T2[n][m][l][k])) << endl;
		for(n=0; n<4; n++) for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k)
			cout << "OMinus_" << "\t" << n+1 << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_OMinus[n][m][l][k])) << "\t" << imag_part((mat_OMinus[n][m][l][k])) << endl;
		for(n=0; n<4; n++) for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k)
			cout << "OPlus_" << "\t" << n+1 << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_OPlus[n][m][l][k])) << "\t" << imag_part((mat_OPlus[n][m][l][k])) << endl;
		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
			cout << "A0_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_A0[m][l][k])) << "\t" << imag_part((mat_A0[m][l][k])) << endl;
		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
			cout << "B0_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_B0[m][l][k])) << "\t" << imag_part((mat_B0[m][l][k])) << endl;
		for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
			cout << "C0_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((mat_C0[m][l][k])) << "\t" << imag_part((mat_C0[m][l][k])) << endl;
*/
}


void fn2F(){
/*
	=---> Screen
*/
int k, l, m, n; // running vars

		cout <<"Quan's Code" << endl << "q_ij" << endl;
/*		cout << evalf(mat_q[0][0]) << "\t" << evalf(mat_q[0][1]) << "\t" << evalf(mat_q[0][2]) << "\t" << evalf(mat_q[0][3]) << "\t" << endl;
		cout << evalf(mat_q[1][0]) << "\t" << evalf(mat_q[1][1]) << "\t" << evalf(mat_q[1][2]) << "\t" << evalf(mat_q[1][3]) << "\t" << endl;
		cout << evalf(mat_q[2][0]) << "\t" << evalf(mat_q[2][1]) << "\t" << evalf(mat_q[2][2]) << "\t" << evalf(mat_q[2][3]) << "\t" << endl;
		cout << evalf(mat_q[3][0]) << "\t" << evalf(mat_q[3][1]) << "\t" << evalf(mat_q[3][2]) << "\t" << evalf(mat_q[3][3]) << "\t" << endl;
*/
		cout << (mat_q[0][0]) << "\t" << (mat_q[0][1]) << "\t" << (mat_q[0][2]) << "\t" << (mat_q[0][3]) << "\t" << endl;
		cout << (mat_q[1][0]) << "\t" << (mat_q[1][1]) << "\t" << (mat_q[1][2]) << "\t" << (mat_q[1][3]) << "\t" << endl;
		cout << (mat_q[2][0]) << "\t" << (mat_q[2][1]) << "\t" << (mat_q[2][2]) << "\t" << (mat_q[2][3]) << "\t" << endl;
		cout << (mat_q[3][0]) << "\t" << (mat_q[3][1]) << "\t" << (mat_q[3][2]) << "\t" << (mat_q[3][3]) << "\t" << endl;
		cout <<"m_i square" << endl;
		cout << (mat_msquare[0]) << "\t" << (mat_msquare[1]) << "\t" << (mat_msquare[2]) << "\t" << (mat_msquare[3]) << "\t" << endl;
for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
	cout << "a_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((fn_a (l, k))) << "\t" << imag_part((fn_a(l, k))) << endl;
for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
	cout << "b_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((fn_b (l, k))) << "\t" << imag_part((fn_b(l, k))) << endl;
for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
	cout << "c_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((fn_c (l, k))) << "\t" << imag_part((fn_c(l, k))) << endl;
for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
	cout << "d_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((fn_d (l, k))) << "\t" << imag_part((fn_d(l, k))) << endl;
for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
	cout << "AC_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((fn_AC (l, k))) << "\t" << imag_part((fn_AC(l, k))) << endl;
for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
	cout << "alpha_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((fn_alpha (l, k))) << "\t" << imag_part((fn_alpha (l, k))) << endl;
for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
	cout << "f_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((fn_f (l, k))) << "\t" << imag_part((fn_f (l, k))) << endl;
for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
	cout << "fminus_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((fn_fminus (l, k))) << "\t" << imag_part((fn_fminus (l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "A_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((fn_A (m, l, k))) << "\t" << imag_part((fn_A (m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "B_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((fn_B (m, l, k))) << "\t" << imag_part((fn_B (m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "C_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((fn_C (m, l, k))) << "\t" << imag_part((fn_C (m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "D_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((fn_D (m, l, k))) << "\t" << imag_part((fn_D (m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "beta_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((fn_beta (m, l, k))) << "\t" << imag_part((fn_beta (m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "phi_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((fn_phi (m, l, k))) << "\t" << imag_part((fn_phi (m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "Q_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((fn_Q (m, l, k))) << "\t" << imag_part((fn_Q (m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "P_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((fn_P (m, l, k))) << "\t" << imag_part((fn_P (m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "E_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((fn_E (m, l, k))) << "\t" << imag_part((fn_E (m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "g_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((fn_g (m, l, k))) << "\t" << imag_part((fn_g (m, l, k))) << endl;	
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "gminus_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((fn_gminus (m, l, k))) << "\t" << imag_part((fn_gminus (m, l, k))) << endl;
for(n=0; n<4; n++) for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k)
	cout << "F_" << "\t" << n+1 << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((fn_F (n, m, l, k))) << "\t" << imag_part((fn_F (n, m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "z1beta_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((fn_z1beta (m, l, k))) << "\t" << imag_part((fn_z1beta (m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "z2beta_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((fn_z2beta (m, l, k))) << "\t" << imag_part((fn_z2beta (m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "z1phi_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((fn_z1phi (m, l, k))) << "\t" << imag_part((fn_z1phi (m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "z2phi_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((fn_z2phi (m, l, k))) << "\t" << imag_part((fn_z2phi (m, l, k))) << endl;
for(n=0; n<4; n++) for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k)
	cout << "T1_" << "\t" << n+1 << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((fn_T1 (n, m, l, k))) << "\t" << imag_part((fn_T1 (n, m, l, k))) << endl;
for(n=0; n<4; n++) for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k)
	cout << "T2_" << "\t" << n+1 << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part((fn_T2 (n, m, l, k))) << "\t" << imag_part((fn_T2 (n, m, l, k))) << endl;

/*
for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
	cout << "a_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_a (l, k))) << "\t" << imag_part(evalf(fn_a(l, k))) << endl;
for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
	cout << "b_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_b (l, k))) << "\t" << imag_part(evalf(fn_b(l, k))) << endl;
for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
	cout << "c_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_c (l, k))) << "\t" << imag_part(evalf(fn_c(l, k))) << endl;
for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
	cout << "d_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_d (l, k))) << "\t" << imag_part(evalf(fn_d(l, k))) << endl;
for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
	cout << "AC_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_AC (l, k))) << "\t" << imag_part(evalf(fn_AC(l, k))) << endl;
for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
	cout << "alpha_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_alpha (l, k))) << "\t" << imag_part(evalf(fn_alpha (l, k))) << endl;
for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
	cout << "f_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_f (l, k))) << "\t" << imag_part(evalf(fn_f (l, k))) << endl;
for(l=0; l<4; l++) for(k=0; k<4; k++) if(l!=k)
	cout << "fminus_" << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_fminus (l, k))) << "\t" << imag_part(evalf(fn_fminus (l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "A_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_A (m, l, k))) << "\t" << imag_part(evalf(fn_A (m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "B_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_B (m, l, k))) << "\t" << imag_part(evalf(fn_B (m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "C_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_C (m, l, k))) << "\t" << imag_part(evalf(fn_C (m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "D_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_D (m, l, k))) << "\t" << imag_part(evalf(fn_D (m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "beta_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_beta (m, l, k))) << "\t" << imag_part(evalf(fn_beta (m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "phi_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_phi (m, l, k))) << "\t" << imag_part(evalf(fn_phi (m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "Q_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_Q (m, l, k))) << "\t" << imag_part(evalf(fn_Q (m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "P_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_P (m, l, k))) << "\t" << imag_part(evalf(fn_P (m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "E_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_E (m, l, k))) << "\t" << imag_part(evalf(fn_E (m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "g_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_g (m, l, k))) << "\t" << imag_part(evalf(fn_g (m, l, k))) << endl;	
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "gminus_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_gminus (m, l, k))) << "\t" << imag_part(evalf(fn_gminus (m, l, k))) << endl;
for(n=0; n<4; n++) for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k)
	cout << "F_" << "\t" << n+1 << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_F (n, m, l, k))) << "\t" << imag_part(evalf(fn_F (n, m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "z1beta_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_z1beta (m, l, k))) << "\t" << imag_part(evalf(fn_z1beta (m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "z2beta_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_z2beta (m, l, k))) << "\t" << imag_part(evalf(fn_z2beta (m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "z1phi_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_z1phi (m, l, k))) << "\t" << imag_part(evalf(fn_z1phi (m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "z2phi_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_z2phi (m, l, k))) << "\t" << imag_part(evalf(fn_z2phi (m, l, k))) << endl;
for(n=0; n<4; n++) for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k)
	cout << "T1_" << "\t" << n+1 << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_T1 (n, m, l, k))) << "\t" << imag_part(evalf(fn_T1 (n, m, l, k))) << endl;
for(n=0; n<4; n++) for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k)
	cout << "T2_" << "\t" << n+1 << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_T2 (n, m, l, k))) << "\t" << imag_part(evalf(fn_T2 (n, m, l, k))) << endl;
*/
cout << evalf(fn_OMinus(0, 1, 2, 3)) << endl;
return;
for(n=0; n<4; n++) for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k)
	cout << "OMinus_" << "\t" << n+1 << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_OMinus (n, m, l, k))) << "\t" << imag_part(evalf(fn_OMinus (n, m, l, k))) << endl;
for(n=0; n<4; n++) for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k)
	cout << "OPlus_" << "\t" << n+1 << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_OPlus (n, m, l, k))) << "\t" << imag_part(evalf(fn_OPlus (n, m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "A0_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_A0 (m, l, k))) << "\t" << imag_part(evalf(fn_A0 (m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "B0_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_B0 (m, l, k))) << "\t" << imag_part(evalf(fn_B0 (m, l, k))) << endl;
for(m=0; m<4; m++) for(l=0; l<4; l++) for(k=0; k<4; k++) if(m!=l && m!=k && l!=k)
	cout << "C0_" << "\t" << m+1 << "\t" << l+1 << "\t" << k+1 << "\t" << real_part(evalf(fn_C0 (m, l, k))) << "\t" << imag_part(evalf(fn_C0 (m, l, k))) << endl;

}
}// Namespace xloops
