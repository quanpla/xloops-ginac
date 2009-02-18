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



	//	III.9	Misc. Functions
	DECLARE_FUNCTION_1P(myfn_delta);//Delta Function
	DECLARE_FUNCTION_1P(myfn_IvZero);
	DECLARE_FUNCTION_1P(myfn_Zero);
	DECLARE_FUNCTION_1P(myfn_IvTheta);
	DECLARE_FUNCTION_1P(myfn_Sign0);
	DECLARE_FUNCTION_2P(fn_eta);

	ex fn_eta_plus(const ex &sigma, const ex &z1, const ex &z2, const ex &P);
	ex fn_eta_minus(const ex &sigma, const ex &z1, const ex &z2, const ex &P);

	ex fn_R(const ex &T1, const ex &T2);
	ex fn_LogACG(const ex &a, const ex &b, const ex &x, const ex &y);
	ex fn_LogARG(const ex &a, const ex &b, const ex &x, const ex &y);
	ex fn_LogAG(const ex &a, const ex &b, const ex &x, const ex &y);

}	// Namespace xloops
#endif 	//  __XLOOPS_ONELOOP_4PT_MY_FUNS_H__
