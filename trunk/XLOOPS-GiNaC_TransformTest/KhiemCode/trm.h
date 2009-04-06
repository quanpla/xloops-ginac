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
**
********************************************************************************
**
** Historial Log:
**	Date		Version	Author		Description
**	_________	_______	_________	________________________________
**	20090216	1.0	Quan Phan	Create this file
**
*******************************************************************************/


#ifndef __XLOOPS_ONELOOP_4PT_TRM_H__
#define __XLOOPS_ONELOOP_4PT_TRM_H__

#include <ginac/ginac.h>
#include <sstream>

namespace xloops{
	using namespace std;
	using namespace GiNaC;
#define EPSILON 0
/*******************************************************************************
	**	I.	Global Variable Declarations
*******************************************************************************/
	int Eps = 0;	// Matrix Dimension = 4 - 2*Eps
	ex input_Rho = 1e-25;	// Feynman prescription
	ex input_Rho1 = 1e-26; // appears in Beta & Phi terms
	ex input_Rho2 = 1e-27; // appears in ThetaG function
	ex Rho = possymbol("Rho"), Rho1 = possymbol("Rho1"), Rho2 = possymbol("Rho2");

	//	I.1	Input Variables
	ex mat_q[4][4], mat_msquare[4];

	//	I.2	Output Variables
	ex D0;	//Final Integrals
/*******************************************************************************
	**	III.	Private Variable Declarations
*******************************************************************************/

	//	III.1	Level 1 Variables
	ex mat_a[4][4], mat_b[4][4];
	ex mat_c[4][4], mat_d[4][4];
	ex mat_d_im[4][4];
	ex mat_d_re[4][4];
	ex mat_d_conj[4][4];

	//	III.2	Level 2 Variables
	ex mat_AC[4][4], mat_alpha[4][4];

	//	III.3	Level 3 Variables
	ex mat_A[4][4][4];
	ex mat_B[4][4][4];
	ex mat_C[4][4][4];
	ex mat_C_im[4][4][4];
	ex mat_C_re[4][4][4];
	ex mat_C_conj[4][4][4];

	ex mat_D[4][4][4];
	ex mat_f[4][4];
	ex mat_fminus[4][4];

	//	III.5	Level 4 Variables
	ex mat_beta[4][4][4];
	ex mat_phi[4][4][4];
	ex mat_g[4][4][4];
	ex mat_gminus[4][4][4];

	//	III.6	Level 5 Variables
	ex mat_Q[4][4][4];
	ex mat_Q_im[4][4][4];
	ex mat_Q_conj[4][4][4];
	ex mat_Q_re[4][4][4];
	ex mat_P[4][4][4];
	ex mat_E[4][4][4];
	ex mat_E_im[4][4][4];
	ex mat_E_re[4][4][4];
	ex mat_F[4][4][4][4];
	ex mat_F_im[4][4][4][4];

	//	III.4	Level 6 Variables
	ex mat_A0[4][4][4];
	ex mat_B0[4][4][4];
	ex mat_C0[4][4][4];
	ex mat_z1phi[4][4][4];
	ex mat_z2phi[4][4][4];
	ex mat_z3phi[4][4][4];
	ex mat_z4phi[4][4][4];
	ex mat_z1beta[4][4][4];
	ex mat_z2beta[4][4][4];
	ex mat_z3beta[4][4][4];
	ex mat_z4beta[4][4][4];

	ex mat_T1[4][4][4][4];
	ex mat_T2[4][4][4][4];
	ex mat_T3[4][4][4][4];
	ex mat_T4[4][4][4][4];

	//	III.7	Level 7 Variables
	ex mat_OPlus[4][4][4][4];
	ex mat_OMinus[4][4][4][4];
	
}	// Namespace xloops

#endif 	// __XLOOPS_ONELOOP_4PT_TRM_H__
