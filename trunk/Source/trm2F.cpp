/** ******************************************************************************
 **      Test Cases
 **1. Write all calculated numeric to file.
 *********************************************************************************/

#include <iostream>
#include <sstream>
#include <fstream>

#include "ginac/ginac.h"
#include "trm2F.h"
#include "extrm.h"

using namespace std;
using namespace GiNaC;

using namespace xloops;

void exportTermsCompareWithKhiem(string filePathName){
 /**
	Export the Result to a file for compare with Khiem code
  */

	ofstream f (filePathName.c_str());

	if (f.is_open()){
		int k, l, m, n;
		f <<"Quan's Code" << endl << "q_ij" << endl;
		f << mat_q[0][0] << "\t" << mat_q[0][1] << "\t" << mat_q[0][2] << "\t" << mat_q[0][3] << "\t" << endl;
		f << mat_q[1][0] << "\t" << mat_q[1][1] << "\t" << mat_q[1][2] << "\t" << mat_q[1][3] << "\t" << endl;
		f << mat_q[2][0] << "\t" << mat_q[2][1] << "\t" << mat_q[2][2] << "\t" << mat_q[2][3] << "\t" << endl;
		f << mat_q[3][0] << "\t" << mat_q[3][1] << "\t" << mat_q[3][2] << "\t" << mat_q[3][3] << "\t" << endl;
		f <<"m_i square" << endl;
		f << mat_msquare[0] << "\t" << mat_msquare[1] << "\t" << mat_msquare[2] << "\t" << mat_msquare[3] << "\t" << endl;
	
		f <<"terms\t real part \t image part" << endl;
		for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) {
			f <<"a_" << k + 1<< l + 1<< "\t" << real_part(mat_a[k][l]) << "\t" << imag_part(mat_a[k][l]) << endl;
		}	
		for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) {
			f <<"b_" << k + 1<< l + 1<< "\t" << real_part(mat_b[k][l]) << "\t" << imag_part(mat_b[k][l]) << endl;
		}		 
		for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) {
			f <<"c_" << k + 1<< l + 1<< "\t" << real_part(mat_c[k][l]) << "\t" << imag_part(mat_c[k][l]) << endl;
		} 
		for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) {
			f <<"d_" << k + 1<< l + 1<< "\t" << real_part(mat_d[k][l]) << "\t" << imag_part(mat_d[k][l]) << endl;
		} 
		for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) {
			f <<"AC_" << k + 1<< l + 1<< "\t" << real_part(mat_AC[k][l]) << "\t" << imag_part(mat_AC[k][l]) << endl;
		} 
		for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) {
			f <<"alpha_" << k + 1<< l + 1<< "\t" << real_part(mat_alpha[k][l]) << "\t" << imag_part(mat_alpha[k][l]) << endl;
		} 
		for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) {
			f <<"f_" << k + 1<< l + 1<< "\t" << real_part(mat_f[k][l]) << "\t" << imag_part(mat_f[k][l]) << endl;
		} 
		for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) {
			f <<"fminus_" << k + 1<< l + 1<< "\t" << real_part(mat_fminus[k][l]) << "\t" << imag_part(mat_fminus[k][l]) << endl;
		} 
		for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) for(m=0; m<4; m++) if(m!=l && m!=k) {
			f <<"A_" << m + 1<< k + 1<< l + 1<< "\t" << real_part(mat_A[m][k][l]) << "\t" << imag_part(mat_A[m][k][l]) << endl;
		} 
		for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) for(m=0; m<4; m++) if(m!=l && m!=k) {
			f <<"B_" << m + 1<< k + 1<< l + 1<< "\t" << real_part(mat_B[m][k][l]) << "\t" << imag_part(mat_B[m][k][l]) << endl;
		} 
		for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) for(m=0; m<4; m++) if(m!=l && m!=k) {
			f <<"C_" << m + 1<< k + 1<< l + 1<< "\t" << real_part(mat_C[m][k][l]) << "\t" << imag_part(mat_C[m][k][l]) << endl;
		} 
		for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) for(m=0; m<4; m++) if(m!=l && m!=k) {
			f <<"D_" << m + 1<< k + 1<< l + 1<< "\t" << real_part(mat_D[m][k][l]) << "\t" << imag_part(mat_D[m][k][l]) << endl;
		} 
		for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) for(m=0; m<4; m++) if(m!=l && m!=k) {
			f <<"beta_" << m + 1<< k + 1<< l + 1<< "\t" << real_part(mat_beta[m][k][l]) << "\t" << imag_part(mat_beta[m][k][l]) << endl;
		} 
		for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) for(m=0; m<4; m++) if(m!=l && m!=k) {
			f <<"phi_" << m + 1<< k + 1<< l + 1<< "\t" << real_part(mat_phi[m][k][l]) << "\t" << imag_part(mat_phi[m][k][l]) << endl;
		} 
		for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) for(m=0; m<4; m++) if(m!=l && m!=k) {
			f <<"Q_" << m + 1<< k + 1<< l + 1<< "\t" << real_part(mat_Q[m][k][l]) << "\t" << imag_part(mat_Q[m][k][l]) << endl;
		}
		for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) for(m=0; m<4; m++) if(m!=l && m!=k) {
			f <<"P_" << m + 1<< k + 1<< l + 1<< "\t" << real_part(mat_P[m][k][l]) << "\t" << imag_part(mat_P[m][k][l]) << endl;
		} 
		for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) for(m=0; m<4; m++) if(m!=l && m!=k) {
			f <<"E_" << m + 1<< k + 1<< l + 1<< "\t" << real_part(mat_E[m][k][l]) << "\t" << imag_part(mat_E[m][k][l]) << endl;
		}
		for(k=0; k<4; k++) for(l=0; l<4; l++)if(l!=k) for(m=0; m<4; m++)if(m!=l && m!=k) {
			f <<"g_" << m + 1<< k + 1<< l + 1<< "\t" << real_part(mat_g[m][k][l]) << "\t" << imag_part(mat_g[m][k][l]) << endl;
		} 
		for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) for(m=0; m<4; m++) if(m!=l && m!=k) {
			f <<"gminus_" << m + 1<< k + 1<< l + 1<< "\t" << real_part(mat_gminus[m][k][l]) << "\t" << imag_part(mat_gminus[m][k][l]) << endl;
		} 
		for(k=0; k<4; k++) for(l=0; l<4; l++) if (l!=k) for(m=0; m<4; m++) if(m!=l && m!=k) for(n=0; n<4; n++) if(n!=l && n!=k && n!=m){
			f <<"F_" << n + 1<< m + 1<< l + 1<< k + 1<< "\t" << real_part(mat_F[n][m][l][k]) << "\t" << imag_part(mat_F[n][m][l][k]) << endl;
		}
	}
	else
		cout << "Cannot write to file !" << endl;
}
