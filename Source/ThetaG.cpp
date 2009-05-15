#include "ginac/ginac.h"
#include "RFunction.h"
#include "ThetaG.h"
#include "trmchk.h" // to raise error when the denominator is zero
using namespace GiNaC;

namespace xloops{
	// calculate Delta, Z1, Z2
	ex fn_delta(const ex &a, const ex &b, const ex&c){
		return pow(b,2) - 4.0*a*c;
	}
	ex fn_z1(const ex &a, const ex &b, const ex &c){
		ex delta = fn_delta(a, b, c);
		return (-b + sqrt(delta))
			/
			(2.0*a);
	}
	ex fn_z2(const ex &a, const ex &b, const ex &c){
		ex delta = fn_delta(a, b, c);
		return (-b - sqrt(delta))
			/
			(2.0*a);
	}

	ex ThetaG(const ex &a, const ex &b, const ex &c, const ex &x, const ex &y){
		ex return_ThetaG;

		// need to check is zero, is positive, is negative, etc.
		int a_pos, a_zero, a_neg, b_pos, b_zero, b_neg, c_pos, c_zero, c_neg, delta_pos, delta_zero, delta_neg;
		a_pos = a.info(info_flags::positive);
		a_zero = a.is_zero();
		a_neg = a.info(info_flags::negative);
		b_pos = b.info(info_flags::positive);
		b_zero = b.is_zero();
		b_neg = b.info(info_flags::negative);
		c_pos = c.info(info_flags::positive);
		c_zero = c.is_zero();
		c_neg = c.info(info_flags::negative);
		delta_pos = fn_delta(a, b, c).info(info_flags::positive);
		delta_zero = fn_delta(a, b, c).is_zero();
		delta_neg = fn_delta(a, b, c).info(info_flags::negative);

		//1.
		if(a_zero && b_zero && !c_neg)
			return_ThetaG = RFunction(-x, -y) + RFunction(x, y);
		//2.
		else if(a_zero && b_zero &&  c_neg)
			return_ThetaG = 0;
		//3.
		else if(a_zero && b_pos)
			return_ThetaG = RFunction (-b/c + x, -b/c + y);
		//4.
		else if(a_zero && b_neg)
			return_ThetaG = RFunction ( b/c - x,  b/c - y);
		//5.
		else if (a_pos && !delta_pos)
			return_ThetaG = RFunction(-x, -y) + RFunction(x, y);
		//6.
		else if (a_pos && delta_pos){
			ex z1 = fn_z1(a, b, c), z2 = fn_z2(a, b, c);
			return_ThetaG = RFunction(-z2 -x, -z2 -y) + RFunction(z1 + x, z1 + y);
		}
		//7.
		else if (a_neg && !delta_pos)
			return_ThetaG = 0;
		//8.
		else if (a_neg && delta_pos){
			ex z1 = fn_z1(a, b, c), z2 = fn_z2(a, b, c);
			return_ThetaG = RFunction(z2 + x, z2 + y) + RFunction(z1 + x, z2 + y);
		}
		return return_ThetaG;
	}
}// Namespace xloops
