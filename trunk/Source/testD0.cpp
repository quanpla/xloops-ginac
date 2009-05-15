#include <iostream>

#include "ginac/ginac.h"
#include "my_fns.h"
#include "D0.h"

using namespace std;
using namespace GiNaC;

int main(){
	ex p, q, m, Rho, Rho1, Rho2;
	ex D0;

	m = lst(6561, 8281, 6561, 8281);
	p = lst(1, 5, 1, 7, 15, 1);
	// convert p --> q
	q = xloops::p2q(p); // definition of the p2q function is in my_fns.h
	Rho = 1e-20; Rho1 = 1e-21; Rho2 = 1e-23;

	ex valueD0 = xloops::D0(q, m, Rho, Rho1, Rho2);
	cout << valueD0.evalf();
	return EXIT_SUCCESS;
}
