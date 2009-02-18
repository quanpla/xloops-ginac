#include <iostream>

#include "ginac/ginac.h"
#include "OneLoop4Pt.h"

using namespace std;
using namespace GiNaC;
using namespace xloops;

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

int main(){
	ex p, q, m, Rho, Rho1, Rho2;
	ex D0;
	
	m = lst(6561, 8281, 6561, 8281);
	p = lst(1, 5, 1, 7, 15, 1);
	// convert p --> q
	q = p2q(p);
	Rho = 1e-20; Rho1 = 1e-21; Rho2 = 1e-23;
	
	cout << "Starting calculating.\n";
	D0 = fn_1Loop4Pt(q, m, Rho, Rho1, Rho2);
	cout << "Finish calculating.\n" << D0.evalf() << endl;
	
	return EXIT_SUCCESS;
}