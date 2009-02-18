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
	 /** *****************************************************************************
		 **	a_lk =     2.(q_l0 - q_k0)
	  *****************************************************n****************************/
		return (mat_q[l][0] - mat_q[k][0])*2.0;
	}

	ex fn_b(int l, int k){
	 /** *****************************************************************************
		 **	b_lk =     -2.(q_l1 - q_k1)
	  *********************************************************************************/
		return -(mat_q[l][1] - mat_q[k][1])*2.0;
	}

	ex fn_c (int l, int k){
	 /** *****************************************************************************
		 **	c_lk =     -2.(q_l2 - q_k2)
	  *********************************************************************************/
		ex c_lk;
		c_lk = -(mat_q[l][2] - mat_q[k][2])*2.0;
		return c_lk;
	}

	ex fn_d (int l, int k){
	 /** *****************************************************************************
		 **	d_lk =     (q_l - q_k)^2  -  (m_l^2 - m_k^2)
	  *********************************************************************************/
		ex d_lk;
		ex temp = 0;
	
		temp = pow(mat_q[l][0] - mat_q[k][0], 2);
	
		for(int ii = 1; ii<4; ii++){
			temp -= pow(mat_q[l][ii]-mat_q[k][ii],2);
		}
	
		d_lk = temp - (mat_msquare[l] - mat_msquare[k]);
		return d_lk;
	}
	ex fn_d_im (int l, int k){
		/** image part of d_lk*/
		return imag_part(mat_d[l][k]);
	}
	ex fn_d_re (int l, int k){
		/** real part of d_lk*/
		return real_part(mat_d[l][k]);
	}
	ex fn_d_conj (int l, int k){
		/** real part of d_lk*/
		return mat_d[l][k].conjugate();
	}

}// Namespace xloops
