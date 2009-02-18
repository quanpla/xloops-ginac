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
**	level 4 functions
**
********************************************************************************
**
** Historial Log:
**	Date		Version	Author		Description
**	_________	_______	_________	________________________________
**	20090214	1.0	Quan Phan	Create this file
**
*******************************************************************************/


#ifndef __XLOOPS_ONELOOP_4PT_LEV4_H__
#define __XLOOPS_ONELOOP_4PT_LEV4_H__

#include <ginac/ginac.h>

using namespace std;

namespace xloops{
	//	III.4	Level 4 Variable Functions
	ex fn_beta (int m, int l, int k);
	ex fn_phi (int m, int l, int k);
	ex fn_g (int m, int l, int k);	
	ex fn_gminus (int m, int l, int k);
	DECLARE_FUNCTION_1P(myfn_g);
	DECLARE_FUNCTION_1P(myfn_gminus);
}	// Namespace xloops

#endif 	// __XLOOPS_ONELOOP_4PT_LEV4_H__
