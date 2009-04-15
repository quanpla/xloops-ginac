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
** Historial Log:
**	Date		Version	Author		Description
**	_________	_______	_________	________________________________
**	20090214	1.0	Quan Phan	Extract this file from the big code
**
*******************************************************************************/


#ifndef __XLOOPS_ONELOOP_4PT_MY_FUNS_H__
#define __XLOOPS_ONELOOP_4PT_MY_FUNS_H__

#include <ginac/ginac.h>

using namespace std;
using namespace GiNaC;

namespace xloops{
	
	/* I override some functions for the case of small argument. Those functions include: is_zero, csgn, step */
	DECLARE_FUNCTION_1P(my_is_zero);
	DECLARE_FUNCTION_1P(my_csgn);
	DECLARE_FUNCTION_1P(my_step);
	DECLARE_FUNCTION_1P(my_is_negative);
	DECLARE_FUNCTION_1P(my_is_positive);

	ex my_evalf(const ex &x);
 
	ex p2q(const ex &p); // convert p to q
}	// Namespace xloops
#endif 	//  __XLOOPS_ONELOOP_4PT_MY_FUNS_H__
