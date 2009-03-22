#include <iostream>
#include "math.h"
#include "ginac/ginac.h"

#include "trm.h"
#include "lev1.h"
#include "lev2.h"
#include "lev3.h"
#include "lev4.h"
#include "logdecmp.h"
// for the function convert p to q
#include "my_fns.h"

using namespace std;
using namespace GiNaC;
using namespace xloops;

void init(const ex& q_, const ex& m_, const ex& Rho_, const ex& Rho1_, const ex& Rho2_){
	ex return_value;

		// init variables
	mat_msquare[0]=m_.op(0); mat_msquare[1]=m_.op(1); mat_msquare[2]=m_.op(2); mat_msquare[3]=m_.op(3);
	mat_q[0][0]=q_.op(0);mat_q[0][1]=0;mat_q[0][2]=0;mat_q[0][3]=0;
	mat_q[1][0]=q_.op(1);mat_q[1][1]=q_.op(2);mat_q[1][2]=0;mat_q[1][3]=0;
	mat_q[2][0]=q_.op(3);mat_q[2][1]=q_.op(4);mat_q[2][2]=q_.op(5);mat_q[2][3]=0;
	input_Rho = Rho_; input_Rho1 = Rho1_; input_Rho2 = Rho2_;
	Rho = Rho_; Rho1 = Rho1_; Rho2 = Rho2_;

		// init terms
	xloopsGiNaC_calc_lev1();
	xloopsGiNaC_calc_lev2();
	xloopsGiNaC_calc_lev3();
	xloopsGiNaC_calc_lev4();
}

int main(int argc, char *argv[]){
	double p10, p20, p21, p30, p31, p32, m1s, m2s, m3s, m4s;
	double xstart = -4000, xstep = 1.0e-10, xend = 0;

	p10 = atof(argv[1]);
	p20 = atof(argv[2]);
	p21 = atof(argv[3]);
	p30 = atof(argv[4]);
	p31 = atof(argv[5]);
	p32 = atof(argv[6]);
	m1s = atof(argv[7]);
	m2s = atof(argv[8]);
	m3s = atof(argv[9]);
	m4s = atof(argv[10]);
	
	ex p, q, m_xloops, Rho, Rho1, Rho2;
	ex theIntegrand;
	
	m_xloops = lst(m1s, m2s, m3s, m4s);
	p = lst(p10, p20, p21, p30, p31, p32);
	// convert p --> q
	q = p2q(p);
	
	Rho = 1e-20; Rho1 = 0; Rho2 = 0;
	
	init(q, m_xloops, Rho, Rho1, Rho2);

	if (argc > 12){
		xstart = atof(argv[11]);
		xstep = atof(argv[12]);
		xend = atof(argv[13]);
	}

	int k, l, m, n;

	// I will loop through the z, check if the value is differ
	double eps = 1e-15;
	// for percent display:
	double one_thoundsandth = fabs(xstart - xend) / 1000.0;
	double current_thoundsandth = 0;

	for(double z = xstart; z <= xend; z += xstep){
		if( z > xstart+current_thoundsandth*one_thoundsandth ){
			while (z > xstart+current_thoundsandth*one_thoundsandth)
				current_thoundsandth += 1;
			printf("\nScanned: %5.1f%%", current_thoundsandth/10.0);
		}
		for (k=0; k<4; k++) for (l=0; l<4; l++) for (m=0; m<4; m++) for (n=0; n<4; n++){
			if(n!=m && n!=l && n!=k && m!=l && m!=k && l!=k){ // compute and compare
				ex log_orig_plus_beta = log_original_plus(k, l, m, n, -1.0/mat_beta[k][l][m], z);
				ex log_orig_minus_beta = log_original_minus(k, l, m, n, -1.0/mat_beta[k][l][m], z);

				ex log_orig_plus_phi = log_original_plus(k, l, m, n, -mat_phi[k][l][m], z);
				ex log_orig_minus_phi = log_original_minus(k, l, m, n, -mat_phi[k][l][m], z);

				ex log_82_plus_beta = log_NPoint_82_plus(k, l, m, n, -1.0/mat_beta[k][l][m], z);
				ex log_82_minus_beta = log_NPoint_82_minus(k, l, m, n, -1.0/mat_beta[k][l][m], z);

				ex log_82_plus_phi = log_NPoint_82_plus(k, l, m, n, -mat_phi[k][l][m], z);
				ex log_82_minus_phi = log_NPoint_82_minus(k, l, m, n, -mat_phi[k][l][m], z);

				ex log_46_plus_beta = log_Khiem_46_plus(k, l, m, n, -1.0/mat_beta[k][l][m], z);
				ex log_46_minus_beta = log_Khiem_46_minus(k, l, m, n, -1.0/mat_beta[k][l][m], z);

				ex log_46_plus_phi = log_Khiem_46_plus(k, l, m, n, -mat_phi[k][l][m], z);
				ex log_46_minus_phi = log_Khiem_46_minus(k, l, m, n, -mat_phi[k][l][m], z);

				if(abs(log_orig_plus_beta - log_46_plus_beta) > eps || abs(log_orig_plus_beta - log_82_plus_beta) > eps){
					printf("\nError at: (k=%d, l=%d, m=%d, n=%d), z = %e\n", k+1, l+1, m+1, n+1, z);
					cout << "z1sigma = " << fn_z1sigma(k, l, m, n, -1.0/mat_beta[k][l][m]) << ", z2sigma = " << fn_z2sigma(k, l, m, n, -1.0/mat_beta[k][l][m]) << endl;
					cout << "Log Orig Plus beta = " << log_orig_plus_beta.evalf() << "\t\tLog 82 Plus beta = " << log_82_plus_beta.evalf() << "\t\tLog 46 Plus beta = " << log_46_plus_beta.evalf() << endl;
				}
				if(abs(log_orig_minus_beta - log_46_minus_beta) > eps || abs(log_orig_minus_beta - log_82_minus_beta) > eps){
					printf("\nError at: (k=%d, l=%d, m=%d, n=%d), z = %e\n", k+1, l+1, m+1, n+1, z);
					cout << "z1sigma = " << fn_z1sigma(k, l, m, n, -1.0/mat_beta[k][l][m]) << ", z2sigma = " << fn_z2sigma(k, l, m, n, -1.0/mat_beta[k][l][m]) << endl;
					cout << "Log Orig Minus beta = " << log_orig_minus_beta.evalf() << "\t\tLog 82 Minus beta = " << log_82_minus_beta.evalf() << "\t\tLog 46 Minus beta = " << log_46_minus_beta.evalf() << endl;
				}

				if(abs(log_orig_plus_phi - log_46_plus_phi) > eps || abs(log_orig_plus_phi - log_82_plus_phi) > eps){
					printf("\nError at: (k=%d, l=%d, m=%d, n=%d), z = %e\n", k+1, l+1, m+1, n+1, z);
					cout << "z1sigma = " << fn_z1sigma(k, l, m, n, -mat_phi[k][l][m]) << ", z2sigma = " << fn_z2sigma(k, l, m, n, -mat_phi[k][l][m]) << endl;
					cout << "Log Orig Plus phi = " << log_orig_plus_phi.evalf() << "\t\tLog 82 Plus phi = " << log_82_plus_phi.evalf() << "\t\tLog 46 Plus phi = " << log_46_plus_phi.evalf() << endl;
				}

				if(abs(log_orig_minus_phi - log_46_minus_phi) > eps || abs(log_orig_minus_phi - log_82_minus_phi) > eps){
					printf("\nError at: (k=%d, l=%d, m=%d, n=%d), z = %e\n", k+1, l+1, m+1, n+1, z);
					cout << "z1sigma = " << fn_z1sigma(k, l, m, n, -mat_phi[k][l][m]) << ", z2sigma = " << fn_z2sigma(k, l, m, n, -mat_phi[k][l][m]) << endl;
					cout << "Log Orig Minus phi = " << log_orig_minus_phi.evalf() << "\t\tLog 82 Minus phi = " << log_82_minus_phi.evalf() << "\t\tLog 46 Minus phi = " << log_46_minus_phi.evalf() << endl;
				}

			}// compute and compare
		} // k, l, m, n run
	} // z runs
	cout << endl;
	return EXIT_SUCCESS;
}
