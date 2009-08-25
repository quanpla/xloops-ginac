#include "ginac/ginac.h"
#include "LogAG.h"
#include "my_fns.h"
#include "trmchk.h" // to raise error when the denominator is zero
using namespace GiNaC;

namespace xloops{
	ex LogAG(const ex &a, const ex &b, const ex &T1, const ex &T2){
		ex rLogAG;
		check0denom(b, "LogAG", 0);
		ex r = a/b;
		rLogAG = ( Li2(1.0 - r*T1) - Li2(1.0 - r*T2) 
					+ myfn_eta(T1, r)*log(1.0 - r*T1) 
					- myfn_eta(T2, r)*log(1.0 - r*T2)
					) / (-T1 + T2)
			+	(log(T1) - log(T2)) * log(b)/(T1-T2);
		return rLogAG;
	}
}// Namespace xloops
