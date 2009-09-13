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
**	level 1 functions
**
********************************************************************************
**
** Historial Log:
**	Date		Version	Author		Description
**	_________	_______	_________	________________________________
**	20090214	1.0	Quan Phan	Create this file
**
*******************************************************************************/


#ifndef __XLOOPS_ONELOOP_4PT_LEV1_H__
#define __XLOOPS_ONELOOP_4PT_LEV1_H__

#include <ginac/ginac.h>

using namespace std;

namespace xloops{
/*******************************************************************************
	**	IV.	Private Function Declarations
*******************************************************************************/
	//	III.1	Level 1 Variable Functions
	ex fn_a (int l, int k);
	ex fn_b (int l, int k);
	ex fn_c (int l, int k);
	ex fn_d (int l, int k); ex fn_d_re(int l, int k); ex fn_d_im(int l, int k); ex fn_d_conj(int l, int k);
	void lev1Calc(); // calculate all the terms (matrix)
	void lev1Calc(int debug); // calculate with debug info
}	// Namespace xloops

#endif 	// __XLOOPS_ONELOOP_4PT_LEV1_H__
