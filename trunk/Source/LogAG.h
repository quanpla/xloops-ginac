#ifndef __XLOOPS_ONELOOP_4PT_LOGAGFUNCTION_H__
#define __XLOOPS_ONELOOP_4PT_LOGAGFUNCTION_H__

#include <ginac/ginac.h>

using namespace std;
using namespace GiNaC;

namespace xloops{
	// LogAG Function from equation [51], Khiem's OneLoop4Pt1.pdf
	ex LogAG(const ex &a, const ex &b, const ex &x, const ex &y);
}	// Namespace xloops
#endif 	//  __XLOOPS_ONELOOP_4PT_LOGAGFUNCTION_H__
