#ifndef __XLOOPS_ONELOOP_4PT_THETAGFUNCTION_H__
#define __XLOOPS_ONELOOP_4PT_THETAGFUNCTION_H__

#include <ginac/ginac.h>

using namespace std;
using namespace GiNaC;

namespace xloops{
	// ThetaGFunction from equation 50, Khiem's OneLoop4Pt1.pdf
	ex ThetaG(const ex &a, const ex &b, const ex &c, const ex &x, const ex &y);
}	// Namespace xloops
#endif 	//  __XLOOPS_ONELOOP_4PT_THETAGFUNCTION_H__
