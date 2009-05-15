/*******************************************************************************
**
**	xloop-ginacs Project
**	HCMUNS,
**	June 2008
**
**	Author(s): 	Son, D. H. 	(sondo01@gmail.com)
**			Khiem, Phan 	(phanhongkhiem@gmail.com)
**			Quan, Phan 	(anhquan.phanle@gmail.com)
********************************************************************************
**
**	Data checking.
**
********************************************************************************
**
**	20090214	1.0	Quan Phan	Extract from the big pile
**
*******************************************************************************/

#ifndef __XLOOPS_ONELOOP_4PT_TRMCHK_H__
#define __XLOOPS_ONELOOP_4PT_TRMCHK_H__

#include <string>
#include "ginac/ginac.h"

using namespace GiNaC;

void check0denom(const ex & denom, std::string varname, int index1, int = -1, int = -1, int = -1);

#endif 	// __XLOOPS_ONELOOP_4PT_TRMCHK_H__
