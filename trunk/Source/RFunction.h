#ifndef __XLOOPS_ONELOOP_4PT_RFUNCTION_H__
#define __XLOOPS_ONELOOP_4PT_RFUNCTION_H__

#include <ginac/ginac.h>

using namespace std;
using namespace GiNaC;

namespace xloops{
	// RFunction from equation 49, Khiem's OneLoop4Pt1.pdf
	ex RFunction(const ex &x, const ex &y);
}	// Namespace xloops
#endif 	//  __XLOOPS_ONELOOP_4PT_RFUNCTION_H__
