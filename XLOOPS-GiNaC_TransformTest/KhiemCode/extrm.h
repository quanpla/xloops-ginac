/*******************************************************************************
**
**	xloop-ginacs Project
**	HCMUNS,
**	June 2008
**
**	Author(s): 	Son, D. H. 	(sondo01@gmail.com)
**			Khiem, Phan 	(phanhongkhiem@gmail.com)
**			Quan, Phan 	(anhquan.phanle@gmail.com)
********************************************************************************
**
**	Extract all terms to file.
**
********************************************************************************
**
**	20090214	1.0	Quan Phan	Extract from the big pile
						This file extern all the terms
**
*******************************************************************************/

#ifndef __XLOOPS_ONELOOP_4PT_EXTRM_H__
#define __XLOOPS_ONELOOP_4PT_EXTRM_H__

#include <ginac/ginac.h>
using namespace GiNaC;

namespace xloops{
/*******************************************************************************
	**	I.	Global Variable Declarations
*******************************************************************************/
	extern int Eps;	// Matrix Dimension = 4 - 2*Eps
	extern ex input_Rho;	// Feynman prescription
	extern ex input_Rho1; // appears in Beta & Phi terms
	extern ex input_Rho2; // appears in ThetaG function
	extern ex Rho, Rho1, Rho2;

	//	I.1	Input Variables
	extern ex mat_q[4][4], mat_msquare[4];



	//	I.2	Output Variables
	extern ex D0;	//Final Integrals


/*******************************************************************************
	**	III.	Private Variable Declarations
*******************************************************************************/

	//	III.1	Level 1 Variables
	extern ex mat_a[4][4], mat_b[4][4];
	extern ex mat_c[4][4], mat_d[4][4];
	extern ex mat_d_im[4][4];
	extern ex mat_d_re[4][4];
	extern ex mat_d_conj[4][4];

	//	III.2	Level 2 Variables
	extern ex mat_AC[4][4], mat_alpha[4][4];

	//	III.3	Level 3 Variables
	extern ex mat_A[4][4][4];
	extern ex mat_B[4][4][4];
	extern ex mat_C[4][4][4];
	extern ex mat_C_im[4][4][4];
	extern ex mat_C_re[4][4][4];
	extern ex mat_C_conj[4][4][4];

	extern ex mat_D[4][4][4];
	extern ex mat_f[4][4];
	extern ex mat_fminus[4][4];

	//	III.5	Level 4 Variables
	extern ex mat_beta[4][4][4];
	extern ex mat_phi[4][4][4];
	extern ex mat_g[4][4][4];
	extern ex mat_gminus[4][4][4];

	//	III.6	Level 5 Variables
	extern ex mat_Q[4][4][4];
	extern ex mat_Q_im[4][4][4];
	extern ex mat_Q_conj[4][4][4];
	extern ex mat_Q_re[4][4][4];
	extern ex mat_P[4][4][4];
	extern ex mat_E[4][4][4];
	extern ex mat_E_im[4][4][4];
	extern ex mat_E_re[4][4][4];
	extern ex mat_F[4][4][4][4];
	extern ex mat_F_im[4][4][4][4];

	//	III.4	Level 6 Variables
	extern ex mat_A0[4][4][4];
	extern ex mat_B0[4][4][4];
	extern ex mat_C0[4][4][4];
	extern ex mat_z1phi[4][4][4];
	extern ex mat_z2phi[4][4][4];
	extern ex mat_z3phi[4][4][4];
	extern ex mat_z4phi[4][4][4];
	extern ex mat_z1beta[4][4][4];
	extern ex mat_z2beta[4][4][4];
	extern ex mat_z3beta[4][4][4];
	extern ex mat_z4beta[4][4][4];

	extern ex mat_T1[4][4][4][4];
	extern ex mat_T2[4][4][4][4];
	extern ex mat_T3[4][4][4][4];
	extern ex mat_T4[4][4][4][4];

	//	III.7	Level 7 Variables
	extern ex mat_OPlus[4][4][4][4];
	extern ex mat_OMinus[4][4][4][4];
}

#endif 	// __XLOOPS_ONELOOP_4PT_EXTRM_H__
