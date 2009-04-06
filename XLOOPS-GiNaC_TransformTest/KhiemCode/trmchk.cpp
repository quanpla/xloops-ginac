#include <iostream>
#include <sstream>
#include <string>

#include "ginac/ginac.h"
#include "trmchk.h"
#include "my_fns.h"

using namespace std;
using namespace GiNaC;

namespace xloops{

	void check0denom(const ex & denom, string varname, int index1, int index2, int index3, int index4){
		// check the input denominator, if it is zero, raise runtime error and equivalent term name
		ex factor = my_evalf(denom);
		if(my_is_zero(factor) == 1.0){
			stringstream err_msg;
			err_msg << varname << "_";
			if (index1 > 0)
				err_msg << index1 + 1;
			if (index2 > 0)
				err_msg << index2 + 1;
			if (index3 > 0)
				err_msg << index3 + 1;
			if (index4 > 0)
				err_msg << index4 + 1;
               
			err_msg << " denominator = 0";
               
			throw std::runtime_error(err_msg.str());
		}
	}
}
