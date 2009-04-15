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
**	Extract all terms to file.
**
********************************************************************************
**
**	20090214	1.0	Quan Phan	Extract from the big pile
						This file extern all the terms
**
*******************************************************************************/

#ifndef __XLOOPS_ONELOOP_4PT_EXTRM_H__
#define __XLOOPS_ONELOOP_4PT_EXTRM_H__

#include <ginac/ginac.h>
using namespace GiNaC;

namespace xloops{
/*******************************************************************************
	**	I.	Global Variable Declarations
*******************************************************************************/
	extern ex Rho, Rho1, Rho2;

	//	I.1	Input Variables
	extern ex mat_q[4][4], mat_msquare[4];

	//	I.2	Output Variables
	extern ex D0;	//Final Integrals
}

#endif 	// __XLOOPS_ONELOOP_4PT_EXTRM_H__
