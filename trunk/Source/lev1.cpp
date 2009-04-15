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
		ex q_l0 = mat_q[l][0], q_k0 = mat_q[k][0];
		
		// calc.
		return 2.0*(q_l0 - q_k0);
	}

	ex fn_b(int l, int k){
		// init
		ex q_l1 = mat_q[l][1], q_k1 = mat_q[k][1];
		
		// calc.
		return 2.0*(q_l1 - q_k1);
	}

	ex fn_c (int l, int k){
		// init
		ex q_l2 = mat_q[l][2], q_k2 = mat_q[k][2];
		
		// calc.
		return 2.0*(q_l2 - q_k2);
	}

	ex fn_d (int l, int k){
		// init
		ex q_l0 = mat_q[l][0], q_l1 = mat_q[l][1], q_l2 = mat_q[l][2], q_l3 = mat_q[l][3],
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
		ex d_lk = fn_d(l, k);
		return imag_part(d_lk);
	}
	ex fn_d_re (int l, int k){
		/** real part of d_lk*/
		ex d_lk = fn_d(l, k);
		return real_part(d_lk);
	}
	ex fn_d_conj (int l, int k){
		/** conjugate of d_lk*/
		ex d_lk = fn_d(l, k);
		return d_lk.conjugate();
	}

	// continue to see lev2.cpp
}// Namespace xloops
