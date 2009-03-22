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

/*20090211: Quan added this function to substitue the Rho and evaluation*/
ex my_evalf(const ex &x){
	ex return_value;

 //      1.      Subs Rho's
	return_value = x.subs(lst(Rho == input_Rho, Rho1 == input_Rho1, Rho2 == input_Rho2));

 //      Last: evaluate
	return_value = return_value.evalf();

	return return_value;
}


void exportTerms2File(string filePathName){
 /**
	Export the Result to a file for Test Case 1
  */

	ofstream f (filePathName.c_str());

	if (f.is_open()){
		int k, l, m, n;
		f << "msquare = ["<< mat_msquare[0] << ","<< mat_msquare[1] << ","<< mat_msquare[2] << ","<< mat_msquare[3] << "]"<<endl;
		f << "q1 = ["<< mat_q[0][0] << ","<< mat_q[0][1] << ","<< mat_q[0][2] << ","<< mat_q[0][3] << "]"<<endl;
		f << "q2 = ["<< mat_q[1][0] << ","<< mat_q[1][1] << ","<< mat_q[1][2] << ","<< mat_q[1][3] << "]"<<endl;
		f << "q3 = ["<< mat_q[2][0] << ","<< mat_q[2][1] << ","<< mat_q[2][2] << ","<< mat_q[2][3] << "]"<<endl;
		f << "rho1 = "<<Rho<<" --- rho2 = "<<Rho2<<endl<<endl;

		f << "----------------------------------------------"<< endl << "1. Level 1 Variables:" << endl;
		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){
			f << "      a_" << k+1 << l+1 << "=";
			if(l!=k) f << mat_a[k][l];
		}
		f << endl;
		}
		f << endl;

		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){
			f << "      b_" << k+1 << l+1 << "=";
			if(l!=k) f << mat_b[k][l];
		}
		f << endl;
		}
		f << endl;

		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){
			f << "      c_" << k+1 << l+1 << "=";
			if(l!=k) f << mat_c[k][l];
		}
		f << endl;
		}
		f << endl;

		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){
			f << "      d_" << k+1 << l+1 << "=";
			if(l!=k) f << mat_d[k][l];
		}
		f << endl;
		}
		f << endl;

		f << "----------------------------------------------"<< endl << "2. Level 2 Variables:" << endl;
		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){
			f << "      AC_" << k+1 << l+1 << "=";
			if(l!=k) f << mat_AC[k][l];
		}
		f << endl;
		}
		f << endl;

		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){
			f << "      alpha_" << k+1 << l+1 << "=";
			if(l!=k) f << mat_alpha[k][l];
		}
		f << endl;
		}
		f << endl;

		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){ for (m = 0; m < 4; m++){
			f << "      A_" << k+1 << l+1 << m+1 << "=";
			if(l!=k && m!=k && l!=m) f << mat_A[k][l][m];
		}
		f << endl;
		}
		f  << endl;
		}
		f << endl;

		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){ for (m = 0; m < 4; m++){
			f << "      B_" << k+1 << l+1 << m+1 << "=";
			if(l!=k && m!=k && l!=m) f << mat_B[k][l][m];
		}
		f << endl;
		}
		f  << endl;
		}
		f << endl;

		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){ for (m = 0; m < 4; m++){
			f << "      C_" << k+1 << l+1 << m+1 << "=";
			if(l!=k && m!=k && l!=m) f << mat_C[k][l][m];
		}
		f << endl;
		}
		f  << endl;
		}
		f << endl;

		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){ for (m = 0; m < 4; m++){
			f << "      D_" << k+1 << l+1 << m+1 << "=";
			if(l!=k && m!=k && l!=m) f << mat_D[k][l][m];
		}
		f << endl;
		}
		f  << endl;
		}
		f << endl;

		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){
			f << "      f_" << k+1 << l+1 << "=";
			if(l!=k) f << mat_f[k][l];
		}
		f << endl;
		}
		f << endl;

		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){
			f << "      fminus_" << k+1 << l+1 << "=";
			if(l!=k) f << mat_fminus[k][l];
		}
		f << endl;
		}
		f << endl;


		f << "----------------------------------------------"<< endl << "4. Level 3 Variables:" << endl;

		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){ for (m = 0; m < 4; m++){
			f << "      beta_" << k+1 << l+1 << m+1 << "=";
			if(l!=k && m!=k && l!=m) f << my_evalf(mat_beta[k][l][m]);
		}
		f << endl;
		}
		f  << endl;
		}
		f << endl;


		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){ for (m = 0; m < 4; m++){
			f << "      phi_" << k+1 << l+1 << m+1 << "=";
			if(l!=k && m!=k && l!=m) f << my_evalf(mat_phi[k][l][m]);
		}
		f << endl;
		}
		f  << endl;
		}
		f << endl;


		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){ for (m = 0; m < 4; m++){
			f << "      g_" << k+1 << l+1 << m+1 << "=";
			if(l!=k && m!=k && l!=m) f << my_evalf(mat_g[k][l][m]);
		}
		f << endl;
		}
		f  << endl;
		}
		f << endl;


		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){ for (m = 0; m < 4; m++){
			f << "      gminus_" << k+1 << l+1 << m+1 << "=";
			if(l!=k && m!=k && l!=m) f << my_evalf(mat_gminus[k][l][m]);
		}
		f << endl;
		}
		f  << endl;
		}
		f << endl;


		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){ for (m = 0; m < 4; m++){
			f << "      Q_" << k+1 << l+1 << m+1 << "=";
			if(l!=k && m!=k && l!=m) f << my_evalf(mat_Q[k][l][m]);
		}
		f << endl;
		}
		f  << endl;
		}
		f << endl;
		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){ for (m = 0; m < 4; m++){
			f << "      Q_im" << k+1 << l+1 << m+1 << "=";
			if(l!=k && m!=k && l!=m) f << my_evalf(mat_Q_im[k][l][m]);
		}
		f << endl;
		}
		f  << endl;
		}
		f << endl;

		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){ for (m = 0; m < 4; m++){
			f << "      Q_re" << k+1 << l+1 << m+1 << "=";
			if(l!=k && m!=k && l!=m) f << my_evalf(mat_Q_re[k][l][m]);
		}
		f << endl;
		}
		f  << endl;
		}
		f << endl;


		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){ for (m = 0; m < 4; m++){
			f << "      P_" << k+1 << l+1 << m+1 << "=";
			if(l!=k && m!=k && l!=m) f << my_evalf(mat_P[k][l][m]);
		}
		f << endl;
		}
		f  << endl;
		}
		f << endl;

		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){ for (m = 0; m < 4; m++){
			f << "      E_" << k+1 << l+1 << m+1 << "=";
			if(l!=k && m!=k && l!=m) f << my_evalf(mat_E[k][l][m]);
		}
		f << endl;
		}
		f  << endl;
		}
		f << endl;

		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){ for (m = 0; m < 4; m++){ for (n = 0; n < 4; n++){
			f << "      F_" << k+1 << l+1 << m+1 << n+1 << "=";
			if(l!=k && m!=k && n!=k && m!=l && n!=l && n!=m) f << my_evalf(mat_F[k][l][m][n]);
		}
		f << endl;
		}
		f  << endl;
		}
		f  << endl;
		}
		f << endl;


		f << "----------------------------------------------"<< endl << "6. Level 4 Variables:" << endl;

		
		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){ for (m = 0; m < 4; m++){
			f << "      z1phi_" << k+1 << l+1 << m+1 << "=";
			if(l!=k && m!=k && l!=m) f << my_evalf(mat_z1phi[k][l][m]);
		}
		f << endl;
		}
		f  << endl;
		}
		f << endl;

		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){ for (m = 0; m < 4; m++){
			f << "      z2phi_" << k+1 << l+1 << m+1 << "=";
			if(l!=k && m!=k && l!=m) f << my_evalf(mat_z2phi[k][l][m]);
		}
		f << endl;
		}
		f  << endl;
		}
		f << endl;

		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){ for (m = 0; m < 4; m++){
			f << "      z3phi_" << k+1 << l+1 << m+1 << "=";
			if(l!=k && m!=k && l!=m) f << my_evalf(mat_z3phi[k][l][m]);
		}
		f << endl;
		}
		f  << endl;
		}
		f << endl;

		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){ for (m = 0; m < 4; m++){
			f << "      z4phi_" << k+1 << l+1 << m+1 << "=";
			if(l!=k && m!=k && l!=m) f << my_evalf(mat_z4phi[k][l][m]);
		}
		f << endl;
		}
		f  << endl;
		}
		f << endl;

		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){ for (m = 0; m < 4; m++){
			f << "      z1beta_" << k+1 << l+1 << m+1 << "=";
			if(l!=k && m!=k && l!=m) f << my_evalf(mat_z1beta[k][l][m]);
		}
		f << endl;
		}
		f  << endl;
		}
		f << endl;

		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){ for (m = 0; m < 4; m++){
			f << "      z2beta_" << k+1 << l+1 << m+1 << "=";
			if(l!=k && m!=k && l!=m) f << my_evalf(mat_z2beta[k][l][m]);
		}
		f << endl;
		}
		f  << endl;
		}
		f << endl;

		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){ for (m = 0; m < 4; m++){
			f << "      z3beta_" << k+1 << l+1 << m+1 << "=";
			if(l!=k && m!=k && l!=m) f << my_evalf(mat_z3beta[k][l][m]);
		}
		f << endl;
		}
		f  << endl;
		}
		f << endl;

		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){ for (m = 0; m < 4; m++){
			f << "      z4beta_" << k+1 << l+1 << m+1 << "=";
			if(l!=k && m!=k && l!=m) f << my_evalf(mat_z4beta[k][l][m]);
		}
		f << endl;
		}
		f  << endl;
		}
		f << endl;

		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){ for (m = 0; m < 4; m++){ for (n = 0; n < 4; n++){
			f << "      T1_" << k+1 << l+1 << m+1 << n+1 << "=";
			if(l!=k && m!=k && n!=k && m!=l && n!=l && n!=m) f << my_evalf(mat_T1[k][l][m][n]);
		}
		f << endl;
		}
		f  << endl;
		}
		f  << endl;
		}
		f << endl;

		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){ for (m = 0; m < 4; m++){ for (n = 0; n < 4; n++){
			f << "      T2_" << k+1 << l+1 << m+1 << n+1 << "=";
			if(l!=k && m!=k && n!=k && m!=l && n!=l && n!=m) f << my_evalf(mat_T2[k][l][m][n]);
		}
		f << endl;
		}
		f  << endl;
		}
		f  << endl;
		}
		f << endl;

		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){ for (m = 0; m < 4; m++){ for (n = 0; n < 4; n++){
			f << "      T3_" << k+1 << l+1 << m+1 << n+1 << "=";
			if(l!=k && m!=k && n!=k && m!=l && n!=l && n!=m) f << my_evalf(mat_T3[k][l][m][n]);
		}
		f << endl;
		}
		f  << endl;
		}
		f  << endl;
		}
		f << endl;

		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){ for (m = 0; m < 4; m++){ for (n = 0; n < 4; n++){
			f << "      T4_" << k+1 << l+1 << m+1 << n+1 << "=";
			if(l!=k && m!=k && n!=k && m!=l && n!=l && n!=m) f << my_evalf(mat_T4[k][l][m][n]);
		}
		f << endl;
		}
		f  << endl;
		}
		f  << endl;
		}
		f << endl;

		f << "----------------------------------------------"<< endl << "7. Level 5 Variables:" << endl;
		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){ for (m = 0; m < 4; m++){ for (n = 0; n < 4; n++){
			f << "      OPlus_" << k+1 << l+1 << m+1 << n+1 << "=";
			if(l!=k && m!=k && n!=k && m!=l && n!=l && n!=m) f << my_evalf(mat_OPlus[k][l][m][n]);
		}
		f << endl;
		}
		f  << endl;
		}
		f  << endl;
		}
		f << endl;


		for(k = 0; k < 4; k++){ for (l = 0; l < 4; l++){ for (m = 0; m < 4; m++){ for (n = 0; n < 4; n++){
			f << "      OMinus_" << k+1 << l+1 << m+1 << n+1 << "=";
			if(l!=k && m!=k && n!=k && m!=l && n!=l && n!=m) f << my_evalf(mat_OMinus[k][l][m][n]);
		}
		f << endl;
		}
		f  << endl;
		}
		f  << endl;
		}
		f << endl;
		f << endl;

		ex tmp = D0.evalf();

		f << "Return:>> D0 = " << tmp << endl;
		f << "divided by I*pi^2 = " << tmp/I/Pi/Pi << endl << endl;

		f.close();
		cout << "Wrote to file!!" << endl;
	}
	else
		cout << "Cannot write to file !" << endl;
}


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
		for(k=0; k<4; k++){
			for(l=0; l<4; l++) if(l!=k) {
				f <<"a_" << k + 1<< l + 1<< "\t" << real_part(mat_a[k][l]) << "\t" << imag_part(mat_a[k][l]) << endl;
			}
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
