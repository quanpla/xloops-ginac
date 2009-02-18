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
**	Level 2 functions
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
#include "lev2.h"

using namespace GiNaC;
namespace xloops{

	//	Level 2 Variable Functions
	ex fn_AC (int l, int k){
	 /**
		 **	AC_lk = a_lk + c_lk
	  **/
		ex AC_lk;
		AC_lk = mat_a[l][k] + mat_c[l][k];
		return AC_lk;
	}

	ex fn_alpha (int l, int k){
		ex alpha_lk;
	 /**
		 **			b_lk
		 **	alpha_lk = -------------
		 **		   a_lk  +  c_lk
		 ** If denominator = 0, we set alpha_lk = 0 by default
	  **/
		ex denom = mat_a[l][k] + mat_c[l][k];
	
		if(my_is_zero(denom) == 1.0){
			alpha_lk = 0.0;
		}
		else{
			alpha_lk = mat_b[l][k] / denom;
		}
		return alpha_lk;
	}

}// Namespace xloops
