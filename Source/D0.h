#ifndef __XLOOPS_ONELOOP_4PT_D0_H__
#define __XLOOPS_ONELOOP_4PT_D0_H__

#include <ginac/ginac.h>

using namespace std;
using namespace GiNaC;

namespace xloops{
	// RFunction from equation 49, Khiem's OneLoop4Pt1.pdf
	ex D0(const ex& q_, const ex& m_, const ex& Rho_, const ex& Rho1_, const ex& Rho2_);
}	// Namespace xloops
#endif 	//  __XLOOPS_ONELOOP_4PT_D0_H__
