#include "ginac/ginac.h"
#include "RFunction.h"
#include "trmchk.h" // to raise error when the denominator is zero
using namespace GiNaC;

namespace xloops{
	ex RFunction(const ex &x, const ex &y){
		check0denom(x-y, "RFunction", 0);
		return 	(log(x) - log(y))
			/
			(x-y);
	}
}// Namespace xloops
