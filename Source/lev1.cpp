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
**	Level 1 functions
**
********************************************************************************
**
** Historial Log:
**	Date		Version	Author		Description
**	_________	_______	_________	________________________________
**	20090214	1.0	Quan Phan	Create this file
*******************************************************************************/

#include "ginac/ginac.h"
#include "extrm.h"
#include "my_fns.h"
#include "lev1.h"

using namespace GiNaC;
namespace xloops{

	//	II.1	Level 1 Variable Functions
	ex fn_a (int l, int k){
		// init
		const ex q_l0 = mat_q[l][0], q_k0 = mat_q[k][0];
		
		// calc.
		return 2.0*(q_l0 - q_k0);
	}

	ex fn_b(int l, int k){
		// init
		const ex q_l1 = mat_q[l][1], q_k1 = mat_q[k][1];
		
		// calc.
		// 20090825, as for OneLoop4PtKhiemMain.cpp use minus sign
		return -2.0*(q_l1 - q_k1);
	}

	ex fn_c (int l, int k){
		// init
		const ex q_l2 = mat_q[l][2], q_k2 = mat_q[k][2];
		
		// calc.
		// 20090825, as for OneLoop4PtKhiemMain.cpp use minus sign
		return -2.0*(q_l2 - q_k2);
	}

	ex fn_d (int l, int k){
		// init
		const ex q_l0 = mat_q[l][0], q_l1 = mat_q[l][1], q_l2 = mat_q[l][2], q_l3 = mat_q[l][3],
			q_k0 = mat_q[k][0], q_k1 = mat_q[k][1], q_k2 = mat_q[k][2], q_k3 = mat_q[k][3],
			msquare_l = mat_msquare[l], msquare_k = mat_msquare[k];

		ex d_lk;
		// calc.
		d_lk = pow(q_l0 - q_k0, 2) - pow(q_l1 - q_k1, 2) - pow(q_l2 - q_k2, 2) - pow(q_l3 - q_k3, 2)
			- (msquare_l - msquare_k);
		
		return d_lk;
	}
	ex fn_d_im (int l, int k){
		/** image part of d_lk*/
		const ex d_lk = fn_d(l, k);
		return imag_part(d_lk);
	}
	ex fn_d_re (int l, int k){
		/** real part of d_lk*/
		const ex d_lk = fn_d(l, k);
		return real_part(d_lk);
	}
	ex fn_d_conj (int l, int k){
		/** conjugate of d_lk*/
		const ex d_lk = fn_d(l, k);
		return d_lk.conjugate();
	}

	void lev1Calc(){
		int l, k;
		for (k = 0; k<4; k++) for (l=0; l<4; l++) if(l!=k){
			mat_a[l][k] = fn_a(l, k);
			mat_b[l][k] = fn_b(l, k);
			mat_c[l][k] = fn_c(l, k);
			mat_d[l][k] = fn_d(l, k);
			mat_d_im[l][k] = fn_d_im(l, k);
			mat_d_re[l][k] = fn_d_re(l, k);
			mat_d_conj[l][k] = fn_d_conj(l, k);
		}
	}
	void lev1Calc(int debug){
		if (debug==0){
			lev1Calc();
			return;
		}
		int l, k;
		printf("Level 1 calculations...\n");
		for (k = 0; k<4; k++) for (l=0; l<4; l++) if(l!=k){
			printf("a[%d,%d]\n", l, k);
			mat_a[l][k] = fn_a(l, k);
			printf("b[%d,%d]\n", l, k);
			mat_b[l][k] = fn_b(l, k);
			printf("c[%d,%d]\n", l, k);
			mat_c[l][k] = fn_c(l, k);
			printf("d[%d,%d]\n", l, k);
			mat_d[l][k] = fn_d(l, k);
			printf("d_im[%d,%d]\n", l, k);
			mat_d_im[l][k] = fn_d_im(l, k);
			printf("d_re[%d,%d]\n", l, k);
			mat_d_re[l][k] = fn_d_re(l, k);
			printf("d_conj[%d,%d]\n", l, k);
			mat_d_conj[l][k] = fn_d_conj(l, k);
		}
	}
	// continue to see lev2.cpp
}// Namespace xloops
