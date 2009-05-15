#include "ginac/ginac.h"
#include "LogAG.h"
#include "my_fns.h"
#include "trmchk.h" // to raise error when the denominator is zero
using namespace GiNaC;

namespace xloops{
	// predefine LogACG, LogARG
	ex fn_A(const ex &a, const ex &b, const ex &x, const ex &y){
		//1. Calc A [54]
		// need to add the requirement: a must be positive
		if (!a.info(info_flags::positive)){
			// raise error here
		}
		ex t = b/a;
		ex A;

		ex Ret = real_part(t), Imt = imag_part(t), Rex = real_part(x), Imx = imag_part(x), Rey = real_part(y), Imy = imag_part(y);

		A = sqrt(pow(Ret, 2) + pow(Imt, 2)) + sqrt(pow(Rex, 2) + pow(Imx, 2)) + sqrt(pow(Rey, 2) + pow(Imy, 2));

		return A;
	}

	ex fn_LogACG(const ex &a, const ex &b, const ex &x, const ex &y){
		ex A = fn_A(a, b, x, y);
		ex x0 = x/A, y0 = y/A, z0 = b/a/A;

		ex LogACG;

		LogACG = (( log(x0)-log(y0) ) / ( A*(x0-y0) )) * log(a*A);
		LogACG += -1.0/( A*(x0-y0) )*(
				 pow(log(x0), 2)/2.0 + pow(log(y0), 2)/2.0 + Li2(1.0-z0/y0) - Li2(1.0-z0/x0)
				+ log(y0) * ( myfn_eta(z0-y0, 1.0/(1.0-y0)) - myfn_eta(z0-y0, -1.0/y0) )
				- log(x0) * ( myfn_eta(z0-x0, 1.0/(1.0-x0)) - myfn_eta(z0-x0, -1.0/x0) )
				+log(1.0-x0/y0)*myfn_eta(z0, 1.0/y0) - log(1.0-z0/x0)*myfn_eta(z0, 1.0/x0)
		);

		return LogACG;
	}

	ex fn_LogARG(const ex &a, const ex &b, const ex &x, const ex &y){
		ex LogARG;

		LogARG = log(b) * (log(x)-log(y)) / (x-y) + fn_LogACG(a/b, 1.0, x, y);
		return LogARG;
	}
	
	ex LogAG(const ex &a, const ex &b, const ex &x, const ex &y){
		ex rLogAG;
		if (a.is_zero()){
			// zero exception
		}
		else if(a.info(info_flags::positive)){
			rLogAG = fn_LogACG(a, b, x, y);
		}
		else if(a.info(info_flags::negative)){
			rLogAG = fn_LogARG(a, b, x, y);
		}
		return rLogAG;
	}
}// Namespace xloops
