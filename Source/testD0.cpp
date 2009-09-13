#include <iostream>

#include "ginac/ginac.h"
#include "my_fns.h"
#include "D0.h"

using namespace std;
using namespace GiNaC;

int main(){
	ex p, q, m, Rho, Rho1, Rho2;
	ex D0;

	m = lst(65610, 82810, 65610, 82810);
	p = lst(10, 50, 10, 70, 150, 10);
	// convert p --> q
	q = xloops::p2q(p); // definition of the p2q function is in my_fns.h
	Rho = 1e-30; Rho1 = 1e-21; Rho2 = 1e-23;
	//cout << q << endl;
	ex valueD0 = xloops::D0(q, m, Rho, Rho1, Rho2);
	cout << valueD0 << endl;
	return EXIT_SUCCESS;
}
