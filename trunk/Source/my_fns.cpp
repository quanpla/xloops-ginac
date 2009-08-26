/*******************************************************************************
**
**	xloop-ginacs Project
**	1Loop4Pt implementation
**
**
**	HCMUNS,
**	June 2008
**
**
**	Author(s): 	Son, D. H. 	(sondo01@gmail.com)
**			Khiem, Phan 	(phanhongkhiem@gmail.com)
**			Quan, Phan 	(anhquan.phanle@gmail.com)
********************************************************************************
**
**	Personal use functions
**
********************************************************************************
**
** Historial Log:
**	Date		Version	Author		Description
**	_________	_______	_________	________________________________
**	20090214	1.0	Quan Phan	Create this file
*******************************************************************************/

#include "ginac/ginac.h"
#include "my_fns.h"
#include "extrm.h"

using namespace GiNaC;

namespace xloops{
	/*20090211: Quan added this function to substitue the Rho and evaluation*/
ex my_evalf(const ex &x){
	ex return_value;

	//	Last: evaluate
	return_value = x.evalf();
	return return_value;
}

ex my_is_zero_eval(const ex &x){
	ex factor = x;
	if(is_a<numeric>(factor)){
		if (abs(factor)<=1e-30)
			return 1.0;
		return 0.0;
	}
	return my_is_zero(x).hold();
}
REGISTER_FUNCTION(my_is_zero, eval_func(my_is_zero_eval));

ex my_csgn_eval(const ex &x){
	ex factor = x;
	if (is_a<numeric>(factor)){
	// for Khiem, x>=0: sign = 1, else -1
		if (factor.info(info_flags::negative) == 1)
			return -1.0;
		return 1.0;
	}
	return my_csgn(x).hold();
}
REGISTER_FUNCTION(my_csgn, eval_func(my_csgn_eval));

ex my_step_eval(const ex &x){
	ex factor = x;

	if (is_a<numeric>(x)){
		if (real_part(x) >= 0.0)
			return 1.0;
		else
			return 0.0;
	}
	return my_step(x).hold();
}
REGISTER_FUNCTION(my_step, eval_func(my_step_eval));


ex my_is_positive_eval(const ex &x){
	ex factor = x;

	if (is_a<numeric>(factor)) {
		if (factor.info(info_flags::positive))
			return 1.0;
		else
			return 0.0;
	}
	return my_step(x).hold();
}
REGISTER_FUNCTION(my_is_positive, eval_func(my_is_positive_eval));

ex my_is_negative_eval(const ex &x){
	ex factor = x;

	if (is_a<numeric>(factor)) {
		if (factor.info(info_flags::negative))
			return 1.0;
		else
			return 0.0;
	}
	return my_step(x).hold();
}
REGISTER_FUNCTION(my_is_negative, eval_func(my_is_negative_eval));

ex myfn_delta_eval(const ex &x){
	if (is_a<numeric>(x)){
		return my_is_zero(x);
	}
	else return myfn_delta(x).hold();
}
REGISTER_FUNCTION(myfn_delta, eval_func(myfn_delta_eval));

ex myfn_eta(const ex &a, const ex &b){
	ex eta;
	ex Ima = imag_part(a), Imb = imag_part(b), Imab = imag_part(a*b), Rea = real_part(a), Reb = real_part(b);
	if (my_is_zero(Ima)==1 && my_is_zero(Imb)==1){
		if (my_is_negative(Rea)==1)
			eta = -2.0 * Pi * I * my_step(Imb);
		else
			eta = 0;
	}
	else if (my_is_zero(Imb)==1 && my_is_zero(Ima)==0){
		if (my_is_negative(Reb)==1)
			eta = -2.0 * Pi * I * my_step(Ima);
		else
			eta = 0;
	}
	else{
		eta = 2.0*Pi*I*( my_step(-Ima)*my_step(-Imb)*my_step(Imab) - my_step(Ima)*my_step(Imb)*my_step(-Imab) );
	}
	return eta;
}

ex p2q(const ex &p){
	ex p12, p23, p13, p1s, p2s, p3s, p4s, p12s, p23s;
	ex q, q10, q20, q21, q30, q31, q32;

	p1s = p.op(0); p2s = p.op(1); p3s = p.op(2); p4s = p.op(3);
	p12s = p.op(4); p23s = p.op(5);

	p12 = (p12s - p1s - p2s)/2.0;
	p23 = (p23s - p2s - p3s)/2.0;
	p13 = (p2s + p4s - p12s - p23s)/2.0;

	q10 = sqrt(p1s);
	q20 = p12/sqrt(p1s) + sqrt(p1s);
	q21 = sqrt( pow(p12/sqrt(p1s) + sqrt(p1s),2) - p12s );
	q30 = p13/sqrt(p1s) + p12/sqrt(p1s) + sqrt(p1s);
	q31 = (q30*q20 - p12s - p13 - p23) / q21;
	q32 = sqrt( pow(q30,2) - pow(q31,2) - p4s );

	q = lst(q10, q20, q21, q30, q31, q32);

	return q;
}

}// Namespace xloops
