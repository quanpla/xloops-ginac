#include <iostream>

#include "ginac/ginac.h"
#include "trm.h"
#include "lev1.h"
#include "lev2.h"
#include "lev3.h"
#include "lev4.h"
#include "lev5.h"
#include "RFunction.h"
#include "ThetaG.h"
#include "LogAG.h"
#include "my_fns.h"

#include "trm2F.h"

using namespace std;
using namespace GiNaC;
using namespace xloops;

void init(const ex& q_, const ex& m_, const ex& Rho_, const ex& Rho1_, const ex& Rho2_){
	mat_msquare[0]=m_.op(0);	mat_msquare[1]=m_.op(1);	mat_msquare[2]=m_.op(2);	mat_msquare[3]=m_.op(3);
	mat_q[0][0]=q_.op(0);mat_q[0][1]=0;mat_q[0][2]=0;mat_q[0][3]=0;
	mat_q[1][0]=q_.op(1);mat_q[1][1]=q_.op(2);mat_q[1][2]=0;mat_q[1][3]=0;
	mat_q[2][0]=q_.op(3);mat_q[2][1]=q_.op(4);mat_q[2][2]=q_.op(5);mat_q[2][3]=0;
	Rho = Rho_; Rho1 = Rho1_; Rho2 = Rho2_;
}

int main(){
	ex p, q, m, Rho, Rho1, Rho2;
	ex D0;
	
	m = lst(65610, 82810, 65610, 82810);
	p = lst(10, 50, 10, 70, 150, 10);
	// convert p --> q
	q = p2q(p); // definition of the p2q function is in my_fns.h
	Rho = 1e-30; Rho1 = 1e-21; Rho2 = 1e-23;
	
	init(q, p, Rho, Rho1, Rho2);
/*	
	// from here you can print any function you want, let have an example:
	for (int l = 0; l<4; l++){for (int k = 0; k<4; k++){
			if (l!=k)
				cout << "a_" << l+1 << k+1 << " = " << fn_a(l, k) << "\t";
			else
				cout << "    #    \t";
		}
		cout << endl;
	}
	for (int l = 0; l<4; l++){for (int k = 0; k<4; k++){
			if (l!=k)
				cout << "b_" << l+1 << k+1 << " = " << fn_b(l, k) << "\t";
			else
				cout << "    #    \t";
		}
		cout << endl;
	}
		for (int l = 0; l<4; l++){for (int k = 0; k<4; k++){
			if (l!=k)
				cout << "c_" << l+1 << k+1 << " = " << fn_c(l, k) << "\t";
			else
				cout << "    #    \t";
		}
		cout << endl;
	}
		for (int l = 0; l<4; l++){for (int k = 0; k<4; k++){
			if (l!=k)
				cout << "d_" << l+1 << k+1 << " = " << fn_d(l, k) << "\t";
			else
				cout << "    #    \t";
		}
		cout << endl;
	}
*/
	cout << LogAG(1.0, 2.0, 4.0, 5.0);
	return EXIT_SUCCESS;
}
