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
**	all the equation to calculate D0 go here
**
********************************************************************************
**
** Historial Log:
**	Date		Version	Author		Description
**	_________	_______	_________	________________________________
**
*******************************************************************************/


#ifndef __XLOOPS_ONELOOP_4PT_D0_H__
#define __XLOOPS_ONELOOP_4PT_D0_H__

#include <ginac/ginac.h>

using namespace GiNaC;

namespace xloops{
	ex D0_integrand(int eqnum, const ex& q_, const ex& m_, const ex& Rho_, const ex& Rho1_, const ex& Rho2_);
}	// Namespace xloops

#endif 	// __XLOOPS_ONELOOP_4PT_D0_H__
