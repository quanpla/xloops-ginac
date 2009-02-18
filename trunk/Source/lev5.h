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
**	level 5 functions
**
********************************************************************************
**
** Historial Log:
**	Date		Version	Author		Description
**	_________	_______	_________	________________________________
**	20090214	1.0	Quan Phan	Create this file
**
*******************************************************************************/


#ifndef __XLOOPS_ONELOOP_4PT_LEV5_H__
#define __XLOOPS_ONELOOP_4PT_LEV5_H__

#include <ginac/ginac.h>

using namespace std;

namespace xloops{

	//	III.5	Level 5 Variable Functions
	ex fn_Q (int m, int l, int k); ex fn_Q_im(int m, int l, int k); ex fn_Q_re(int m, int l, int k); ex fn_Q_conj(int m, int l, int k);
	ex fn_P (int m, int l, int k);
	ex fn_E (int m, int l, int k); ex fn_E_im(int m, int l, int k); ex fn_E_re(int m, int l, int k);
	ex fn_F(int n, int m, int l, int k); ex fn_F_im(int n, int m, int l, int k);

}	// Namespace xloops

#endif 	// __XLOOPS_ONELOOP_4PT_LEV5_H__
