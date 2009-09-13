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
	ex fn_z1phi (int m, int l, int k);
	ex fn_z2phi (int m, int l, int k);
	ex fn_z1beta (int m, int l, int k);
	ex fn_z2beta (int m, int l, int k);

	ex fn_T1 (int n, int m, int l, int k);
	ex fn_T2 (int n, int m, int l, int k);

	void lev4Calc();
	void lev4Calc(int debug);	
}	// Namespace xloops

#endif 	// __XLOOPS_ONELOOP_4PT_LEV4_H__
