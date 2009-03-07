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
**	OneLoop4Pt main entry functions
**
********************************************************************************
**
** Historial Log:
**	Date		Version	Author		Description
**	_________	_______	_________	________________________________
**	20090214	1.0	Quan Phan	Create this file
**
*******************************************************************************/


#ifndef __XLOOPS_ONELOOP_4PT_H__
#define __XLOOPS_ONELOOP_4PT_H__

#include <ginac/ginac.h>

using namespace std;
using namespace GiNaC;

namespace xloops{
	//	III.8	Level 8 Function
	ex fn_1Loop4Pt();
	ex fn_1Loop4Pt(const ex& q_, const ex& m_, const ex& Rho, const ex& Rho1, const ex& Rho2);
}	// Namespace xloops

#endif 	// __XLOOPS_ONELOOP_4PT_H__