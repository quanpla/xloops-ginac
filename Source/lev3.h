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
**	level 3 functions
**
********************************************************************************
**
** Historial Log:
**	Date		Version	Author		Description
**	_________	_______	_________	________________________________
**	20090214	1.0	Quan Phan	Create this file
**
*******************************************************************************/


#ifndef __XLOOPS_ONELOOP_4PT_LEV3_H__
#define __XLOOPS_ONELOOP_4PT_LEV3_H__

#include <ginac/ginac.h>

using namespace std;

namespace xloops{

	//	III.3	Level 3 Variable Functions
	ex fn_A (int m, int l, int k);
	ex fn_B (int m, int l, int k);
	ex fn_C (int m, int l, int k); ex fn_C_im (int m, int l, int k); ex fn_C_re (int m, int l, int k); ex fn_C_conj (int m, int l, int k);
	ex fn_D (int m, int l, int k);
	ex fn_f (int l, int k);
	ex fn_fminus (int l, int k);
	DECLARE_FUNCTION_1P(myfn_f);//because there is if condition
	DECLARE_FUNCTION_1P(myfn_fminus);//because there is if condition

}	// Namespace xloops

#endif 	// __XLOOPS_ONELOOP_4PT_LEV3_H__
