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
**
********************************************************************************
**
** Historial Log:
**	Date		Version	Author		Description
**	_________	_______	_________	________________________________
**	20090216	1.0	Quan Phan	Create this file
**
*******************************************************************************/


#ifndef __XLOOPS_ONELOOP_4PT_TRM_H__
#define __XLOOPS_ONELOOP_4PT_TRM_H__

#include <ginac/ginac.h>
#include <sstream>

namespace xloops{
	using namespace std;
	using namespace GiNaC;
#define EPSILON 0
/*******************************************************************************
	**	I.	Global Variable Declarations
*******************************************************************************/

	ex Rho, Rho1, Rho2;

	//	I.1	Input Variables
	ex mat_q[4][4], mat_msquare[4];


}	// Namespace xloops

#endif 	// __XLOOPS_ONELOOP_4PT_TRM_H__
