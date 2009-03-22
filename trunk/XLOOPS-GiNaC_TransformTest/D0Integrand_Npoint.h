/*******************************************************************************
**
**      xloop-ginacs Project
**      For testing purpose, check NPoints equations
**
**
**      HCMUNS,
**      Mar 200
**
**
**      Author(s):      Son, D. H.      (sondo01@gmail.com)
**                      Khiem, Phan     (phanhongkhiem@gmail.com)
**                      Quan, Phan      (anhquan.phanle@gmail.com)
********************************************************************************
**
** Historial Log:
**      Date            Version Author          Description
**      _________       _______ _________       ________________________________
*******************************************************************************/

#ifndef __XLOOPS_ONELOOP_4PT_NPOINTTEST_H__
#define __XLOOPS_ONELOOP_4PT_NPOINTTEST_H__

#include <ginac/ginac.h>

using namespace GiNaC;

namespace xloops{
	ex NPoint_Equation(int eqnum, const ex& q_, const ex& m_, const ex& Rho_, const ex& Rho1_, const ex& Rho2_);
}	// Namespace xloops

#endif 	// __XLOOPS_ONELOOP_4PT_NPOINTTEST_H__

