#include <iostream>

#include "ginac/ginac.h"

#include "ginac/ginac.h"
#include "D0Integrand.h"
// for the function convert p to q
#include "my_fns.h"

using namespace GiNaC;
using namespace xloops;

int main(int argc, char *argv[]){
	int equationNumber = 0;
	double p10, p20, p21, p30, p31, p32, m1s, m2s, m3s, m4s;
	
	equationNumber = atoi(argv[1]);
	p10 = atof(argv[2]);
	p20 = atof(argv[3]);
	p21 = atof(argv[4]);
	p30 = atof(argv[5]);
	p31 = atof(argv[6]);
	p32 = atof(argv[7]);
	m1s = atof(argv[8]);
	m2s = atof(argv[9]);
	m3s = atof(argv[10]);
	m4s = atof(argv[11]);
	
	ex p, q, m, Rho, Rho1, Rho2;
	ex theIntegrand;
	
	m = lst(m1s, m2s, m3s, m4s);
	p = lst(p10, p20, p21, p30, p31, p32);
	// convert p --> q
	q = p2q(p);
	
	Rho = 1e-20; Rho1 = 1e-21; Rho2 = 1e-23;
	
	theIntegrand = D0_integrand(equationNumber, q, m, Rho, Rho1, Rho2);
		
	cout << csrc << real_part(theIntegrand);
	
	return EXIT_SUCCESS;
}