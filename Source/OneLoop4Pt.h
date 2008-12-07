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
**	Notes:
**	Build command: g++ OneLoop4Pt.cpp -lginac -o OneLoop4Pt.out
**	Please see the adjoint documents for detail descriptions
**
********************************************************************************
**
** Historial Log:
**	Date		Version	Author		Description
**	_________	_______	_________	________________________________
**	20080804	1.0	Quan Phan	Add imag_part for d, C, E, Q, Q*, F,
**						z1beta, z2beta, z1phi, z2phi
**	20080715	1.0	Quan Phan	Add fn_eta_beta and fn_eta_phi
**	20080715	1.0	Quan Phan	Create the interface for CMD Line
**	20080531	1.0	Quan Phan	Create this file
**
*******************************************************************************/

/*******************************************************************************
 **  xloops Copyright (C) 1997,2000-2002 Johannes Gutenberg University Mainz,
 **  Germany
 **
 **  This program is free software; you can redistribute it and/or modify
 **  it under the terms of the GNU General Public License as published by
 **  the Free Software Foundation; either version 2 of the License, or
 **  (at your option) any later version.
 **
 **  This program is distributed in the hope that it will be useful,
 **  but WITHOUT ANY WARRANTY; without even the implied warranty of
 **  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 **  GNU General Public License for more details.
 **
 **  You should have received a copy of the GNU General Public License
 **  along with this program; if not, write to the Free Software
 **  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *******************************************************************************/


#ifndef __XLOOPS_ONELOOP_4PT_H__
#define __XLOOPS_ONELOOP_4PT_H__

#include <ginac/ginac.h>
#include <sstream>

//namespace xloops{
	using namespace std;
	using GiNaC::ex;
	using GiNaC::lst;
	using GiNaC::numeric;
	using GiNaC::symbol;
	#define EPSILON 0
/*******************************************************************************
**	I.	Global Variable Declarations
*******************************************************************************/
	int Eps = 0;	// Matrix Dimension = 4 - 2*Eps
	ex Rho = 1e-25;	// Feynman prescription
	ex Rho1 = 1e-26; // appears in Beta & Phi terms
	ex Rho2 = 1e-27; // appears in ThetaG function

	//	I.1	Input Variables
	ex mat_q[4][4], mat_msquare[4];
	
	
	
	//	I.2	Output Variables
	ex D0;	//Final Integrals

	
	
/*******************************************************************************
**	II.	Global Function Declarations
*******************************************************************************/
	
	// The ultimate function, implementation of the Execution Plan.
//ex fn_1Loop4Pt(lst P, lst M, int Eps, float Rho);
	ex fn_1Loop4Pt();
	ex fn_1Loop4Pt(const ex& q_, const ex& m_, const ex& Rho, const ex& Rho1, const ex& Rho2);
	
	
/*******************************************************************************
**	III.	Private Variable Declarations
*******************************************************************************/

	//	III.1	Level 1 Variables
	ex mat_a[4/*-2*EPSILON*/][4/*-2*EPSILON*/], mat_b[4/*-2*EPSILON*/][4/*-2*EPSILON*/];
	ex mat_c[4/*-2*EPSILON*/][4/*-2*EPSILON*/], mat_d[4/*-2*EPSILON*/][4/*-2*EPSILON*/];
		ex mat_d_im[4/*-2*EPSILON*/][4/*-2*EPSILON*/];
		ex mat_d_re[4/*-2*EPSILON*/][4/*-2*EPSILON*/];
		ex mat_d_conj[4/*-2*EPSILON*/][4/*-2*EPSILON*/];
	
	//	III.2	Level 2 Variables
	ex mat_AC[4/*-2*EPSILON*/][4/*-2*EPSILON*/], mat_alpha[4/*-2*EPSILON*/][4/*-2*EPSILON*/];
	
	//	III.3	Level 3 Variables
	ex mat_A[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
	ex mat_B[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
	ex mat_C[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
		ex mat_C_im[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
		ex mat_C_re[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
		ex mat_C_conj[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];

	ex mat_D[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
	ex mat_f[4/*-2*EPSILON*/][4/*-2*EPSILON*/];
	ex mat_fminus[4/*-2*EPSILON*/][4/*-2*EPSILON*/];
	
	//	III.5	Level 4 Variables
	ex mat_beta[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
	ex mat_phi[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
	ex mat_g[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
	ex mat_gminus[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
	
	//	III.6	Level 5 Variables
	ex mat_Q[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
		ex mat_Q_im[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
		ex mat_Q_conj[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
		ex mat_Q_re[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
	ex mat_P[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
	ex mat_E[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
		ex mat_E_im[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
		ex mat_E_re[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
	ex mat_F[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
		ex mat_F_im[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];

	//	III.4	Level 6 Variables
	ex mat_A0[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
	ex mat_B0[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
	ex mat_C0[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
	ex mat_z1phi[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
		ex mat_z1phi_im[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
	ex mat_z2phi[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
		ex mat_z2phi_im[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
	ex mat_z1beta[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
		ex mat_z1beta_im[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
	ex mat_z2beta[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
		ex mat_z2beta_im[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];

	ex mat_T1[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
	ex mat_T2[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
	
	//	III.7	Level 7 Variables
	ex mat_OPlus[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
	ex mat_OMinus[4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/][4/*-2*EPSILON*/];
	
	
	
/*******************************************************************************
**	IV.	Private Function Declarations
*******************************************************************************/
	//	III.1	Level 1 Variable Functions
	ex fn_a (int l, int k);
	ex fn_b (int l, int k);
	ex fn_c (int l, int k);
	ex fn_d (int l, int k); ex fn_d_im(int l, int k);
	
	//	III.2	Level 2 Variable Functions
	ex fn_AC (int l, int k);
	ex fn_alpha (int l, int k);
	
	//	III.3	Level 3 Variable Functions
	ex fn_A (int m, int l, int k);
	ex fn_B (int m, int l, int k);
	ex fn_C (int m, int l, int k);ex fn_C_im (int m, int l, int k);
	ex fn_D (int m, int l, int k);
	ex fn_f (int l, int k);
	ex fn_fminus (int l, int k);
	DECLARE_FUNCTION_1P(myfn_f);//because there is if condition
	DECLARE_FUNCTION_1P(myfn_fminus);//because there is if condition
	
	//	III.4	Level 4 Variable Functions
	ex fn_beta (int m, int l, int k);
	ex fn_phi (int m, int l, int k);
	ex fn_g (int m, int l, int k);	
	ex fn_gminus (int m, int l, int k);
	DECLARE_FUNCTION_1P(myfn_g);
	DECLARE_FUNCTION_1P(myfn_gminus);
	
	//	III.5	Level 5 Variable Functions
	ex fn_Q (int m, int l, int k); ex fn_Q_im(int m, int l, int k); ex fn_Q_re(int m, int l, int k); ex fn_Q_conj(int m, int l, int k);
	ex fn_P (int m, int l, int k);
	ex fn_E (int m, int l, int k); ex fn_E_im(int m, int l, int k); ex fn_E_re(int m, int l, int k);
	ex fn_F(int n, int m, int l, int k); ex fn_F_im(int n, int m, int l, int k);
	
	//	III.6	Level 6 Variable Functions
	ex fn_A0 (int m, int l, int k);
	ex fn_B0 (int m, int l, int k);
	ex fn_C0 (int m, int l, int k);
	ex fn_z1phi (int m, int l, int k); ex fn_z1phi_im (int m, int l, int k);
	ex fn_z2phi (int m, int l, int k); ex fn_z2phi_im (int m, int l, int k);
	ex fn_z1beta (int m, int l, int k); ex fn_z1beta_im (int m, int l, int k);
	ex fn_z2beta (int m, int l, int k); ex fn_z2beta_im (int m, int l, int k);
	ex fn_T1 (int n, int m, int l, int k);
	ex fn_T2 (int n, int m, int l, int k);
	
	//	III.7	Level 7 Variable Functions
	ex fn_OPlus (int n, int m, int l, int k);
	ex fn_OMinus (int n, int m, int l, int k);
	
	//	III.8	Level 8 Function
	void fn_1Loop4Pt_calc();// Use to calculate 1Loop4Pt after having inputed

	//	III.9	Misc. Functions
	DECLARE_FUNCTION_1P(myfn_delta);//Delta Function
	DECLARE_FUNCTION_1P(myfn_IvZero);
	DECLARE_FUNCTION_1P(myfn_Zero);
	DECLARE_FUNCTION_1P(myfn_IvTheta);
	DECLARE_FUNCTION_1P(myfn_Sign0);

	DECLARE_FUNCTION_2P(fn_eta);

	ex fn_R(const ex &T1, const ex &T2);
//	ex fn_eta(const ex &a, const ex &b);
	ex fn_eta_beta(int m, int l, int k);
	ex fn_eta_phi(int m, int l, int k);
	ex fn_GPositive(const ex &T1, const ex &T2);
	ex fn_GNegative(const ex &T1, const ex &T2);
	ex fn_LnGPositive(const ex &T1, const ex &T2, const ex &z0);
	ex fn_LnGNegative(const ex &T1, const ex &T2, const ex &z0);
	ex fn_ThetaG(const ex &A0, const ex &B0, const ex &C0, const ex &T1, const ex &T2);
//}	// Namespace xloops
	
	//	IV.1. TestCases
	void exportTerms2File(string filePathName);
#endif 	// __XLOOPS_ONELOOP_4PT_H__
