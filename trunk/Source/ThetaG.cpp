#include "ginac/ginac.h"
#include "RFunction.h"
#include "ThetaG.h"
#include "trmchk.h" // to raise error when the denominator is zero
using namespace GiNaC;

namespace xloops{
	ex ThetaG(const ex &a, const ex &b, const ex &x, const ex &y){
		ex return_ThetaG;

		// need to check is zero, is positive, is negative, etc.
		int a_pos, a_zero, a_neg, b_pos, b_zero, b_neg;
		a_pos = a.info(info_flags::positive);
		a_zero = a.is_zero();
		a_neg = a.info(info_flags::negative);
		b_pos = b.info(info_flags::positive);
		b_zero = b.is_zero();
		b_neg = b.info(info_flags::negative);

		//1.
		if (a_zero) {
			if (!b_neg){ //  b >= 0
				return_ThetaG = RFunction(x, y) + RFunction(-x, -y);
			}
			else { // b < 0
				return_ThetaG = 0;
			}
		}
		else if (a_pos){
			return_ThetaG = RFunction(-b/a + x, -b/a + y);
		}
		else{ // a_neg
			return_ThetaG = RFunction(b/a - x, b/a - y);
		}
		return return_ThetaG;
	}

	ex ThetaGC(const ex &a, const ex &b, const ex &c, const ex &x, const ex &y){
		ex return_ThetaGC;

		// need to check is zero, is positive, is negative, etc.
		ex delta = b*b - 4.0*a*c;
		
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
		delta_pos = delta.info(info_flags::positive);
		delta_zero = delta.is_zero();
		delta_neg = delta.info(info_flags::negative);

		//1.
		if (a_zero && b_zero){
			if (!c_pos){
				return_ThetaGC = RFunction(-x, -y) + RFunction(x, y);
			}
			else{
				return_ThetaGC = 0;
			}
		}
		else if (b_pos)
			return_ThetaGC = RFunction(-c/b + x, -c/b + y);
		else if (b_neg)
			return_ThetaGC = RFunction(c/b - x, c/b - y);
		else if (a_pos && !delta_pos)
			return_ThetaGC = RFunction(x, y) + RFunction(-x, -y);
		else if (a_pos && delta_pos)
			return_ThetaGC = 0;
		else if (a_neg && delta_pos){
			ex 	z1 = (-b + sqrt(delta)) / (2.0*a),
				z2 = (-b - sqrt(delta)) / (2.0*a);
			return_ThetaGC = RFunction(z2 + x, z2 + y) + RFunction(z1 + x, z1 + y);
		}
		return return_ThetaGC;
	}
}// Namespace xloops
