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
**	Notes:
**	Build command:
**	 g++ `pkg-config --cflags --libs ginac` OneLoop4Pt.h OneLoop4Pt.cpp -o run.exe
**
**	Please see the adjoint documents for detail descriptions
**	Sucessfully test with ginac 1.3.7
**
********************************************************************************
**
** Historial Log:
**	Date		Version	Author		Description
**	_________	_______	_________	________________________________
**	20081105	1.0	Quan Phan	Clone this file
*******************************************************************************/

 /********************************************************************************
 **
 **	Execution Plan
 **			+-------------------------------+
 **			|				|
 **			| Input P, Q, epsilon, rho 	|
 **			|				|
 **			+-------------------------------+
 **					|
 **					V
 **			+-------------------------------+
 **			|				|
 **			| Symbolic Calculate Level 1	|
 **			| Intermidate Terms 	|
 **			|				|
 **			+-------------------------------+
 **   |
 **   V
 **   ...
 **
 **   |
 **   V
 **			+-------------------------------+
 **			|				|
 **			| Symbolic Calculate Level 7	|
 **			| Intermidate Terms 	|
 **			|				|
 **			+-------------------------------+
 **   |
 **   V
 **			+-------------------------------+
 **			|				|
 **			| Symbolic Calculate Final Terms|
 **			|				|
 **			+-------------------------------+
 **   |
 **   V
 **			+-------------------------------+
 **			|				|
 **			| Subtitute numerical		|
 **			|	 Final Term		|
 **			|				|
 **			+-------------------------------+
 **
 *******************************************************************************/

#include <iostream>
#include <sstream>
#include <fstream>
#include <getopt.h>

#include "ginac/ginac.h"
#include "ginac/flags.h"
#include "OneLoop4Pt_QuanCodeLastEquation.h"

/*include for Vegas*/
/**/
#include <complex>
extern "C" {
#include "vegas.h"
}

/*defines for Vegas*/
// Number of integral dimensions
#define DIMENSION 3

// The integrand need a reference variable, hence this dimension
#define FUNCTIONS 1



using namespace std;
using namespace GiNaC;


/*******************************************************************************
**	I.	Demo Execution Start Point
*******************************************************************************/
void input_manual(lst &q, lst &m, ex &Rho, ex &Rho1, ex &Rho2);
void input_option(int argc, char *argv[], lst &q, lst &m, ex &Rho, ex &Rho1, ex &Rho2);
void output_CompareWithKhiem();
void calculate_D0_Integrand();
void getD0Integrand_Vegas(double x[DIMENSION], double f[FUNCTIONS]);
void output_Subtitue1P();

// global for Vegas to calculate
ex 	tan_piX = symbol("tan_piX");/*tan (pi*x/2)*/
ex	tan_piX_s = symbol("tan_piX_s"); /*tan^2(pi*x/2)*/

static int imagCompute_flag; // =0: compute real part, = 1: compute imaginary part


int main(int argc, char *argv[]){
	lst m, q;
	ex Rho = 0, Rho1=0, Rho2=0;
	ex result;
	if(argc > 1)
	{
		input_option(argc, argv, q, m, Rho, Rho1, Rho2);
		result = fn_1Loop4Pt(q, m, Rho, Rho1, Rho2);
	}
	else
	{
		char question_;
		cout << "Do you want to input NUMERIC values for m & q (y/n):";
		cin >> question_;
		if(question_=='y' || question_=='Y'){
			input_manual(q, m, Rho, Rho1, Rho2);
			result = fn_1Loop4Pt(q, m, Rho, Rho1, Rho2);
		}
		else{
			result = fn_1Loop4Pt();
		}
	}
	//output_CompareWithKhiem();
	calculate_D0_Integrand();
	
	// print out the values of the integrand
	output_Subtitue1P();
	
	return EXIT_SUCCESS;
}

void calculate_D0_Integrand(){
	//the loop will be 4 indices
	int k, l, m, n;
	// init
	D0 = 0;
	
	for(k=0; k<4; k++) for(l=0; l<4; l++) if (l!=k) for(m=0; m<4; m++) if(m!=l && m!=k) for(n=0; n<4; n++) if(n!=l && n!=k && n!=m){
		/*equation 23, luxury computations*/
		ex G_tan = 1.0/(
		                 mat_beta[m][l][k] * ( mat_E[m][l][k]*tan_piX - mat_msquare[k] + I*Rho) /*Feynman's prescription*/
		                 -
		                 (mat_Q[m][l][k] + mat_P[m][l][k]*tan_piX) * (mat_F[n][m][l][k] + tan_piX)
		               );
		ex G_tan_minus = 1.0/(
		                       mat_beta[m][l][k] * (-mat_E[m][l][k]*tan_piX - mat_msquare[k] + I*Rho)/*Feynman's Prescription*/
		                       -
		                       (mat_Q[m][l][k] - mat_P[m][l][k]*tan_piX) * (mat_F[n][m][l][k] - tan_piX)
		                     );
		/*equation 23, 1st line*/
		ex outside_constant;
		outside_constant = 
			mat_B[m][l][k]/(
			                 mat_B[m][l][k]*mat_A[n][l][k]
			                 -
			                 mat_B[n][l][k]*mat_A[m][l][k]
			               )
			*
			(1.0 - myfn_delta(mat_AC[l][k]))
			*
			(1.0 - myfn_delta(mat_B[m][l][k]))
			*
			abs(1.0 - mat_beta[m][l][k]*mat_phi[m][l][k]);
		
		/*equ 23, 2nd to 5th line*/ 
		ex first_integrand;
		first_integrand = 
			(1.0 + tan_piX_s) * G_tan
			*(
			   2.0*log( mat_F[n][m][l][k]/mat_beta[m][l][k] )
			   -
			   log(
			        (
			          (1.0 - mat_beta[m][l][k]*mat_phi[m][l][k]) * tan_piX
			          +
			          mat_F[n][m][l][k]
			        )
			        /
			        mat_beta[m][l][k]
			      )
			   /*line 3*/
			   -
			   log(
			        -(
			           (1.0 - mat_beta[m][l][k]*mat_phi[m][l][k])*tan_piX
			           + mat_F[n][m][l][k]
			         )
			        / mat_beta[m][l][k]
			      )
			   /*line 4*/
			   -
			   2.0 * log(
			              (
			                - mat_P[m][l][k]/mat_beta[m][l][k]*tan_piX_s
			                + (mat_E[m][l][k] - mat_Q[m][l][k]/mat_beta[m][l][k]) * tan_piX
			                - mat_msquare[k] + I*Rho
			              )
			              /(mat_Q[m][l][k] + mat_P[m][l][k]*tan_piX)
			            )
			   /*line 5*/
			   + log(
			          (
			            - mat_P[m][l][k]*mat_phi[m][l][k]*tan_piX_s
			            + (mat_E[m][l][k] - mat_Q[m][l][k]/mat_beta[m][l][k])*tan_piX
			            - mat_msquare[k] + I*Rho /*Feynman*/
			          )
			          /
			         (mat_Q[m][l][k] + mat_P[m][l][k]*tan_piX)
			        )
			   /*line 6*/
			   + log(
			          (
			            -mat_P[m][l][k]*mat_phi[m][l][k]*tan_piX_s
			            +(mat_E[m][l][k] - mat_Q[m][l][k]*mat_phi[m][l][k])*tan_piX
			            -mat_msquare[k] + I*Rho
			          )
			         /(mat_Q[m][l][k] + mat_P[m][l][k]*tan_piX)
			        )
			 );
		/*equ 23, 7th to 11th line*/ 
		ex second_integrand;
		second_integrand = 
			(1.0 + tan_piX_s) * G_tan_minus
			*(
			   - log(mat_F[n][m][l][k]/mat_beta[m][l][k]) - log(-mat_F[n][m][l][k]/mat_beta[m][l][k])
			   /*line 8*/
			   + 2.0*log(
			              (
			                -(1.0 - mat_beta[m][l][k]*mat_phi[m][l][k])*tan_piX + mat_F[n][m][l][k]
			              )
			              /
			              mat_beta[m][l][k]
			            )
			   /*line 9*/
			   + log(
			          (
			            -mat_P[m][l][k]/mat_beta[m][l][k]*tan_piX_s
			            -(mat_E[m][l][k] - mat_Q[m][l][k]/mat_beta[m][l][k])*tan_piX
			            -mat_msquare[k] + I*Rho
			          )
			          /(mat_Q[m][l][k] - mat_P[m][l][k]*tan_piX)
			        )
			   /*line 10*/
			   + log(
			          (
			            -mat_P[m][l][k]/mat_beta[m][l][k]*tan_piX_s
			            -(mat_E[m][l][k] - mat_Q[m][l][k]/mat_beta[m][l][k])*tan_piX
			            -mat_msquare[k] + I*Rho
			          )
			          /(mat_Q[m][l][k] - mat_P[m][l][k]*tan_piX)
			        )
			   /*line 11*/
			   - 2.0 * log(
			                (
			                  - mat_P[m][l][k]*mat_phi[m][l][k]*tan_piX_s
			                  - (mat_E[m][l][k] - mat_Q[m][l][k]*mat_phi[m][l][k])*tan_piX
			                  - mat_msquare[k] + I*Rho /*Feynman*/
			                )
			                /
			                ( mat_Q[m][l][k] -  mat_P[m][l][k]*tan_piX)
			              )
			 );
		D0 += outside_constant * (first_integrand + second_integrand);
	}
	D0 = D0*pow(Pi,3)*I/2.0;
	D0 = D0.evalf();
}

void getD0Integrand_Vegas(double x[DIMENSION], double f[FUNCTIONS]){
	ex out_D0 = 0.0;
	//* substitutions variables
	double 	s_tan_piX = ex_to<numeric>(tan(Pi*x[0]/2.0).evalf()).to_double(),
		s_tan_piX_s = ex_to<numeric>(pow(tan(Pi*x[0]/2.0), 2).evalf()).to_double();
	
	out_D0 = D0.subs(lst(tan_piX==s_tan_piX, tan_piX_s==s_tan_piX_s));
	out_D0 = out_D0.evalf();
	
	if (is_a<numeric>(out_D0)){
		// check if the calculate is a numeric, else raise error
		if (imagCompute_flag == 0)
			f[0] = ex_to<numeric>(real_part(out_D0)).to_double();
		else
			f[0] = ex_to<numeric>(imag_part(out_D0)).to_double();
	}
	else{
		stringstream err_msg;
		err_msg << "Cannot evaluation D0 to float. D0 = " << out_D0 << endl;
		throw runtime_error(err_msg.str());
	}
}

void output_Subtitue1P(){
	for (double x0 = 0.0; x0 <= 1.0; x0 += 0.02){
		double x[DIMENSION] = {x0},
			f[FUNCTIONS];
		imagCompute_flag = 0; // compute the real part
		getD0Integrand_Vegas(x,f);
		cout << x0 << "\t" << f[0];
		imagCompute_flag = 1;
		getD0Integrand_Vegas(x,f);
		cout << "\t" << f[0] << endl;
	}
}

void output_CompareWithKhiem(){
	int k, l, m, n;
	
	cout << "Quan's Code" << endl << "q_ij" << endl;
	for(k=0; k<4; k++){
		for(l=0; l<4; l++)
			cout << mat_q[k][l] << "\t";
		cout << endl;
	}
	
	cout << "m_i square" << endl;
	for(k=0; k<4; k++){
		cout << mat_msquare[k] << "\t";
	}
	cout << endl;
	
	cout << "terms\treal part\timage part" << endl;
	for(k=0; k<4; k++){
		for(l=0; l<4; l++) if(l!=k) {
			cout << "a_" << k+1 << l+1 << "\t" << real_part(mat_a[k][l]) << "\t" << imag_part(mat_a[k][l]) << endl;
		}
	}	
	for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) {
		cout << "b_" << k+1 << l+1 << "\t" << real_part(mat_b[k][l]) << "\t" << imag_part(mat_b[k][l]) << endl;
	}		 
	for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) {
		cout << "c_" << k+1 << l+1 << "\t" << real_part(mat_c[k][l]) << "\t" << imag_part(mat_c[k][l]) << endl;
	} 
	for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) {
		cout << "d_" << k+1 << l+1 << "\t" << real_part(mat_d[k][l]) << "\t" << imag_part(mat_d[k][l]) << endl;
	} 
	for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) {
		cout << "AC_" << k+1 << l+1 << "\t" << real_part(mat_AC[k][l]) << "\t" << imag_part(mat_AC[k][l]) << endl;
	} 
	for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) {
		cout << "alpha_" << k+1 << l+1 << "\t" << real_part(mat_alpha[k][l]) << "\t" << imag_part(mat_alpha[k][l]) << endl;
	} 
	for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) {
		cout << "f_" << k+1 << l+1 << "\t" << real_part(mat_f[k][l]) << "\t" << imag_part(mat_f[k][l]) << endl;
	} 
	for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) {
		cout << "fminus_" << k+1 << l+1 << "\t" << real_part(mat_fminus[k][l]) << "\t" << imag_part(mat_fminus[k][l]) << endl;
	} 
	for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) for(m=0; m<4; m++) if(m!=l && m!=k) {
		cout << "A_" << m+1 << k+1 << l+1 << "\t" << real_part(mat_A[m][k][l]) << "\t" << imag_part(mat_A[m][k][l]) << endl;
	} 
	for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) for(m=0; m<4; m++) if(m!=l && m!=k) {
		cout << "B_" << m+1 << k+1 << l+1 << "\t" << real_part(mat_B[m][k][l]) << "\t" << imag_part(mat_B[m][k][l]) << endl;
	} 
	for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) for(m=0; m<4; m++) if(m!=l && m!=k) {
		cout << "C_" << m+1 << k+1 << l+1 << "\t" << real_part(mat_C[m][k][l]) << "\t" << imag_part(mat_C[m][k][l]) << endl;
	} 
	for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) for(m=0; m<4; m++) if(m!=l && m!=k) {
		cout << "D_" << m+1 << k+1 << l+1 << "\t" << real_part(mat_D[m][k][l]) << "\t" << imag_part(mat_D[m][k][l]) << endl;
	} 
	for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) for(m=0; m<4; m++) if(m!=l && m!=k) {
		cout << "beta_" << m+1 << k+1 << l+1 << "\t" << real_part(mat_beta[m][k][l]) << "\t" << imag_part(mat_beta[m][k][l]) << endl;
	} 
	for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) for(m=0; m<4; m++) if(m!=l && m!=k) {
		cout << "phi_" << m+1 << k+1 << l+1 << "\t" << real_part(mat_phi[m][k][l]) << "\t" << imag_part(mat_phi[m][k][l]) << endl;
	} 
	for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) for(m=0; m<4; m++) if(m!=l && m!=k) {
		cout << "Q_" << m+1 << k+1 << l+1 << "\t" << real_part(mat_Q[m][k][l]) << "\t" << imag_part(mat_Q[m][k][l]) << endl;
	}
	for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) for(m=0; m<4; m++) if(m!=l && m!=k) {
		cout << "P_" << m+1 << k+1 << l+1 << "\t" << real_part(mat_P[m][k][l]) << "\t" << imag_part(mat_P[m][k][l]) << endl;
	} 
	for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) for(m=0; m<4; m++) if(m!=l && m!=k) {
		cout << "E_" << m+1 << k+1 << l+1 << "\t" << real_part(mat_E[m][k][l]) << "\t" << imag_part(mat_E[m][k][l]) << endl;
	}
	for(k=0; k<4; k++) for(l=0; l<4; l++) for(m=0; m<4; m++) {
		cout << "g_" << m+1 << k+1 << l+1 << "\t" << real_part(mat_g[m][k][l]) << "\t" << imag_part(mat_g[m][k][l]) << endl;
	} 
	for(k=0; k<4; k++) for(l=0; l<4; l++) if(l!=k) for(m=0; m<4; m++) if(m!=l && m!=k) {
		cout << "gminus_" << m+1 << k+1 << l+1 << "\t" << real_part(mat_gminus[m][k][l]) << "\t" << imag_part(mat_gminus[m][k][l]) << endl;
	} 
	for(k=0; k<4; k++) for(l=0; l<4; l++) if (l!=k) for(m=0; m<4; m++) if(m!=l && m!=k) for(n=0; n<4; n++) if(n!=l && n!=k && n!=m){
		cout << "F" << n+1 << m+1 << l+1 << k+1 << "\t" << real_part(mat_F[n][m][l][k]) << "\t" << imag_part(mat_F[n][m][l][k]) << endl;
	}
}

//namespace xloops{
ex fn_1Loop4Pt(const ex& q_, const ex& m_, const ex& Rho_, const ex& Rho1_, const ex& Rho2_){
	mat_msquare[0]=m_.op(0);	mat_msquare[1]=m_.op(1);	mat_msquare[2]=m_.op(2);	mat_msquare[3]=m_.op(3);
	mat_q[0][0]=q_.op(0);mat_q[0][1]=0;mat_q[0][2]=0;mat_q[0][3]=0;
	mat_q[1][0]=q_.op(1);mat_q[1][1]=q_.op(2);mat_q[1][2]=0;mat_q[1][3]=0;
	mat_q[2][0]=q_.op(3);mat_q[2][1]=q_.op(4);mat_q[2][2]=q_.op(5);mat_q[2][3]=0;
	Rho = Rho_; Rho1 = Rho1_; Rho2 = Rho2_;
	fn_1Loop4Pt_calc();
	return D0;
}

ex fn_1Loop4Pt(){
	mat_msquare[0]=symbol("m1");mat_msquare[1]=symbol("m2");mat_msquare[2]=symbol("m3");mat_msquare[3]=symbol("m4");
	mat_q[0][0]=symbol("q10");mat_q[0][1]=0;mat_q[0][2]=0;mat_q[0][3]=0;
	mat_q[1][0]=symbol("q20");mat_q[1][1]=symbol("q21");mat_q[1][2]=0;mat_q[1][3]=0;
	mat_q[2][0]=symbol("q30");mat_q[2][1]=symbol("q31");mat_q[2][2]=symbol("q32");mat_q[2][3]=0;
	
	fn_1Loop4Pt_calc();
	return D0;
}

/* I override some functions for the case of small argument. Those functions include: is_zero, csgn, step */
DECLARE_FUNCTION_1P(my_is_zero);
DECLARE_FUNCTION_1P(my_csgn);
DECLARE_FUNCTION_1P(my_step);

ex my_is_zero_eval(const ex &x){
	ex factor = x;
	if(is_a<numeric>(factor)){
		if (abs(factor)<=1e-10)
			return 1;
		return 0;
	}
	return my_is_zero(x).hold();
}
REGISTER_FUNCTION(my_is_zero, eval_func(my_is_zero_eval));

ex my_csgn_eval(const ex &x){
	ex factor = x;
	if (is_a<numeric>(factor)){
		if (my_is_zero(factor) == 1)
			return 0;
		return csgn(x);
	}
	return my_csgn(x).hold();
}
REGISTER_FUNCTION(my_csgn, eval_func(my_csgn_eval));

ex my_step_eval(const ex &x){
	if (is_a<numeric>(x)){
		if (real_part(x) >= 0)
			return 1;
		else
			return 0;
	}
	return my_step(x).hold();
}
REGISTER_FUNCTION(my_step, eval_func(my_step_eval));


	/*******************************************************************************
	**	II.	Functions Implementation
	*******************************************************************************/
	//	II.1	Level 1 Variable Functions
ex fn_a (int l, int k){
	 /** *****************************************************************************
	 **	a_lk = 2.(q_l0 - q_k0)
	 *****************************************************n****************************/
	return (mat_q[l][0] - mat_q[k][0])*2;
}

ex fn_b(int l, int k){
	 /** *****************************************************************************
	 **	b_lk = -2.(q_l1 - q_k1)
	 *********************************************************************************/
	return -(mat_q[l][1] - mat_q[k][1])*2;
}

ex fn_c (int l, int k){
	 /** *****************************************************************************
	 **	c_lk = -2.(q_l2 - q_k2)
	 *********************************************************************************/
	ex c_lk;
	c_lk = -(mat_q[l][2] - mat_q[k][2])*2;
	return c_lk;
}

ex fn_d (int l, int k){
	 /** *****************************************************************************
	 **	d_lk = (q_l - q_k)^2 - (m_l^2 - m_k^2)
	 *********************************************************************************/
	ex d_lk;
	ex temp = 0;
	
	temp = pow(mat_q[l][0] - mat_q[k][0], 2);
	
	for(int ii = 1; ii<4; ii++){
		temp -= pow(mat_q[l][ii]-mat_q[k][ii],2);
	}
	
	d_lk = temp - (mat_msquare[l] - mat_msquare[k]);
	return d_lk;
}
ex fn_d_im (int l, int k){
	 /** image part of d_lk*/
	return imag_part(mat_d[l][k]);
}
ex fn_d_re (int l, int k){
		 /** real part of d_lk*/
	return real_part(mat_d[l][k]);
}
ex fn_d_conj (int l, int k){
		 /** real part of d_lk*/
	return mat_d[l][k].conjugate();
}

	//	II.2	Level 2 Variable Functions
ex fn_AC (int l, int k){
	 /**
	 **	AC_lk = a_lk + c_lk
	 **/
	ex AC_lk;
	AC_lk = mat_a[l][k] + mat_c[l][k];
	return AC_lk;
}

ex fn_alpha (int l, int k){
	ex alpha_lk;
	 /**
	 **			b_lk
	 **	alpha_lk = -------------
	 **		 a_lk + c_lk
	 ** If denominator = 0, we set alpha_lk = 0 by default
	 **/
	ex denom = mat_a[l][k] + mat_c[l][k];
	
	if(my_is_zero(denom)==1){
		alpha_lk = 0;
	}
	else{
		alpha_lk = mat_b[l][k] / denom;
	}
	return alpha_lk;
}

	//	II.3	Level 3 Variable Functions
ex fn_A (int m, int l, int k){
	ex A_mlk;
	 /**
	 **			 AC_mk
	 **	A_mlk = a_mk - a_lk -------
	 **			 AC_lk
	 **/
	ex denom = mat_AC[l][k];
	if(my_is_zero(denom)==1){
		stringstream err_msg;
		err_msg << "A_" << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	
	A_mlk = mat_a[m][k] - mat_a[l][k]*mat_AC[m][k]/mat_AC[l][k];
	return A_mlk;
}

ex fn_B (int m, int l, int k){
	ex B_mlk;
	 /**
	 **			 AC_mk
	 **	B_mlk = b_mk - b_lk -------
	 **			 AC_lk
	 **/
	ex denom = mat_AC[l][k];
	if(my_is_zero(denom)==1){
		stringstream err_msg;
		err_msg << "B_" << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	
	B_mlk = mat_b[m][k] - mat_b[l][k]*mat_AC[m][k]/mat_AC[l][k];
	return B_mlk;
}

ex fn_C (int m, int l, int k){
	ex C_mlk;
	 /**
	 **			 AC_mk
	 **	C_mlk = d_mk - d_lk -------
	 **			 AC_lk
	 **/
	ex denom = mat_AC[l][k];
	if(my_is_zero(denom)==1){
		stringstream err_msg;
		err_msg << "C_" << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	
	C_mlk = mat_d[m][k] - mat_d[l][k]*mat_AC[m][k]/mat_AC[l][k];
	return C_mlk;
}
ex fn_C_im (int m, int l, int k){
	ex C_mlk_im;
	 /**
	 **			 AC_mk
	 **	C_mlk_im = d_mk_im - d_lk_im -------
	 **			 AC_lk
	 **/
	if (is_a<numeric>(mat_C[m][l][k]))
		C_mlk_im = imag_part(mat_C[m][l][k]);
	else{
		ex denom = mat_AC[l][k];
		if(my_is_zero(denom)==1){
			stringstream err_msg;
			err_msg << "C_" << m+1 << l+1 << k+1 << " denominator = 0";
			throw std::runtime_error(err_msg.str());
		}
		C_mlk_im = mat_d_im[m][k] - mat_d_im[l][k]*mat_AC[m][k]/denom;
	}
	return C_mlk_im;
}
ex fn_C_re (int m, int l, int k){
	ex C_mlk_re;
	 /**
	 **			 AC_mk
	 **	C_mlk_re = d_mk_re - d_lk_re -------
	 **			 AC_lk
	 **/
	if (is_a<numeric>(mat_C[m][l][k]))
		C_mlk_re = real_part(mat_C[m][l][k]);
	else{
		ex denom = mat_AC[l][k];
		if(my_is_zero(denom)==1){
			stringstream err_msg;
			err_msg << "C_" << m+1 << l+1 << k+1 << " denominator = 0";
			throw std::runtime_error(err_msg.str());
		}
		C_mlk_re = mat_d_re[m][k] - mat_d_re[l][k]*mat_AC[m][k]/denom;
	}
	return C_mlk_re;
}
ex fn_C_conj (int m, int l, int k){
	ex C_mlk_conj;
	 /**
	 **			  AC_mk
	 **	C_mlk_conj = d_mk_conj - d_lk_conj -------
	 **			  AC_lk
	 **/
	if (is_a<numeric>(mat_C[m][l][k]))
		C_mlk_conj = mat_C[m][l][k].conjugate();
	else{
		ex denom = mat_AC[l][k];
		if(my_is_zero(denom)==1){
			stringstream err_msg;
			err_msg << "C_" << m+1 << l+1 << k+1 << " denominator = 0";
			throw std::runtime_error(err_msg.str());
		}
		C_mlk_conj = mat_d_conj[m][k] - mat_d_conj[l][k]*mat_AC[m][k]/denom;
	}
	return C_mlk_conj;
}

ex fn_D (int m, int l, int k){
	ex D_mlk;
	 /**
	 **		 (q_l - q_k)^2
	 **	D_mlk = - 4 ---------------
	 **			AC_lk^2
	 **/
	ex denom = mat_AC[l][k]*mat_AC[l][k];
	if(my_is_zero(denom)==1){
		stringstream err_msg;
		err_msg << "D_" << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	ex temp = 0;
	temp = pow(mat_q[l][0] - mat_q[k][0], 2);
	for(int ii = 1; ii < 4; ii++){
		ex temp2 = mat_q[l][ii] - mat_q[k][ii];
		temp -= temp2*temp2;
	}
	D_mlk = -(temp / denom) * 4;
	return D_mlk;
}


ex myfn_f_eval(const ex &factor){
	if(is_a<numeric>(factor)){
		if (my_is_zero(factor)==1)
			return 1;
		if (factor.info(info_flags::negative))
			return 0;
		return 2;
	}
	else{
		return myfn_f(factor).hold();
	}
}
REGISTER_FUNCTION(myfn_f, eval_func(myfn_f_eval));

ex fn_f (int l, int k){
	 /**		/-	 d_lk
	 **	 	| = 0, if Im( - ------- ) < 0;
	 **		|		 AC_lk
	 **		|
	 **	 	|	 d_lk
	 **	f_lk	{ = 1, if Im( - ------- ) = 0;
	 **	 	|		 AC_lk
	 **		|
	 **		|	 d_lk
	 **		| = 2, if Im( - ------- ) > 0;
	 **		|		 AC_lk
	 **		\-
	 **/
	ex denom = mat_AC[l][k];
	if(my_is_zero(denom)==1){
		stringstream err_msg;
		err_msg << "f_" << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	
	ex factor = - mat_d_im[l][k]/denom;
	return myfn_f(factor);
}

ex myfn_fminus_eval(const ex &factor){
	if(is_a<numeric>(factor)){
		if (my_is_zero(factor)==1)
			return 1;
		if (factor.info(info_flags::negative))
			return 2;
		return 0;
	}
	else{
		return myfn_fminus(factor).hold();
	}
}

REGISTER_FUNCTION(myfn_fminus, eval_func(myfn_fminus_eval));
ex fn_fminus (int l, int k){
	 /**			/-	 d_lk
	 **	 		| = 0, if Im( - ------- ) > 0;
	 **			|		 AC_lk
	 **			|
	 **	 	 	|	 d_lk
	 **	fminus_lk	{ = 1, if Im( - ------- ) = 0;
	 **	 		|		 AC_lk
	 **			|
	 **			|	 d_lk
	 **			| = 2, if Im( - ------- ) < 0;
	 **			|		 AC_lk
	 **			\-
	 **/
	ex denom = mat_AC[l][k];
	
	if(my_is_zero(denom)==1){
		stringstream err_msg;
		err_msg << "fminus_" << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	ex factor = - mat_d_im[l][k]/denom;
	return myfn_fminus(factor);
}

	//	II.4	Level 4 Variable Functions
ex fn_beta (int m, int l, int k){
	 /**	 A_mlk
	 ** term1 = ------- - alpha_lk
	 **	 B_mlk
	 **
	 **		term1 + sqrt(term1^2 - D_mlk + i*rho)
	 ** beta_mlk = ----------------------------------------
	 **				D_mlk
	 **/
	ex beta_mlk;
	ex denom = mat_D[m][l][k];
	if(my_is_zero(mat_D[m][l][k])==1|| my_is_zero(mat_B[m][l][k])==1){
		stringstream err_msg;
		err_msg << "beta_" << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	
	ex term1 = mat_A[m][l][k] / mat_B[m][l][k] - mat_alpha[l][k];
	
	beta_mlk = term1 + sqrt(term1*term1 - mat_D[m][l][k] + I * Rho1 /*This is NOT Feynman's prescription.*/);
	beta_mlk = beta_mlk / mat_D[m][l][k];
	
	return beta_mlk;
}

ex fn_phi (int m, int l, int k){
	 /**	 A_mlk
	 ** term1 = ------- - alpha_lk
	 **	 B_mlk
	 **
	 ** beta_mlk = term1 + sqrt(term1^2 - D_mlk + i*eta)
	 **/
	ex phi_mlk;
	ex denom = mat_B[m][l][k];
	if (my_is_zero(denom)==1){
		stringstream err_msg;
		err_msg << "phi_" << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	
	ex term1 = mat_A[m][l][k] / mat_B[m][l][k] - mat_alpha[l][k];
	
	phi_mlk = term1 + sqrt(term1*term1 - mat_D[m][l][k] + I * Rho1 /*This is NOT Feynman's prescription*/);
	
	return phi_mlk;
}

ex myfn_g_eval(const ex &factor){
	if (is_a<numeric>(factor)){
		if (my_is_zero(factor)==1)
			return 1;
		if (factor.info(info_flags::negative))
			return 0;
		return 2;
	}
	else{
		return myfn_g(factor).hold();
	}
}

REGISTER_FUNCTION(myfn_g, eval_func(myfn_g_eval));
ex fn_g (int m, int l, int k){
	 /**			/-	 -C_mlk
	 **	 		| = 0, if Im( -------- ) < 0;
	 **			|		B_mlk
	 **			|
	 **	 	 	|	 -C_mlk
	 **	g_mlk		| = 1, if Im( -------- ) = 0;
	 **	 		|		B_mlk
	 **			|
	 **			|	 -C_mlk
	 **			| = 2, if Im( -------- ) > 0;
	 **			|		B_mlk
	 **			\-
	 **/
	ex denom = mat_B[m][l][k];
	if(my_is_zero(denom)==1){
		stringstream err_msg;
		err_msg << "g_" << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	ex factor = -mat_C_im[m][l][k]/denom;
	return myfn_g(factor);
}

ex myfn_gminus_eval(const ex &factor){
	if (is_a<numeric>(factor)){
		if (my_is_zero(factor)==1)
			return 1;
		if (factor.info(info_flags::negative))
			return 2;
		return 0;
	}
	else{
		return myfn_gminus(factor).hold();
	}
}

REGISTER_FUNCTION(myfn_gminus, eval_func(myfn_gminus_eval));
ex fn_gminus (int m, int l, int k){
	 /**			/-	 -C_mlk
	 **	 		| = 0, if Im( -------- ) > 0;
	 **			|		B_mlk
	 **			|
	 **	 	 	|	 -C_mlk
	 **	gminus_mlk	{ = 1, if Im( -------- ) = 0;
	 **	 		|		B_mlk
	 **			|
	 **			|	 -C_mlk
	 **			| = 2, if Im( -------- ) < 0;
	 **			|		B_mlk
	 **			\-
	 **/
	ex denom = mat_B[m][l][k];
	if(my_is_zero(denom)==1){
		stringstream err_msg;
		err_msg << "gminus_" << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	
	ex factor = -mat_C_im[m][l][k]/denom;
	
	return myfn_gminus(factor);
}

	//	II.5	Level 5 Variable Functions
ex fn_Q (int m, int l, int k){
	ex Q_mlk;
	 /** C_mlk d_lk
	 ** Q_mlk = -2 ( ------- + ------- beta_mlk )
	 ** B_mlk AC_lk
	 **/
	if (my_is_zero(mat_B[m][l][k])==1 || my_is_zero(mat_AC[l][k])==1){
		stringstream err_msg;
		err_msg << "Q_" << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	Q_mlk = -(mat_C[m][l][k] / mat_B[m][l][k] + mat_beta[m][l][k] * mat_d[l][k] / mat_AC[l][k]) * 2;
	return Q_mlk;
}
ex fn_Q_im (int m, int l, int k){
	ex Q_mlk_im;
	 /**  Im(C_mlk) Im(d_lk)
	 ** Im(Q_mlk) = -2 ( ---------- + ---------- beta_mlk )
	 **  B_mlk AC_lk
	 **/
	if (my_is_zero(mat_B[m][l][k])==1 || my_is_zero(mat_AC[l][k])==1){
		stringstream err_msg;
		err_msg << "Q_" << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	
	if (is_a<numeric>(mat_Q[m][l][k]))
		Q_mlk_im = imag_part(mat_Q[m][l][k]);
	else{
		Q_mlk_im = mat_C_im[m][l][k]/mat_B[m][l][k] + mat_d_im[l][k]*mat_beta[m][l][k]/mat_AC[l][k];
		Q_mlk_im *= -2;
	}
	return Q_mlk_im;
}
ex fn_Q_re (int m, int l, int k){
	ex Q_mlk_re;
	 /**  Re(C_mlk) Re(d_lk)
	 ** Re(Q_mlk) = -2 ( ---------- + ---------- beta_mlk )
	 **  B_mlk AC_lk
	 **/
	if (my_is_zero(mat_B[m][l][k])==1 || my_is_zero(mat_AC[l][k])==1){
		stringstream err_msg;
		err_msg << "Q_" << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	if (is_a<numeric>(mat_Q[m][l][k]))
		Q_mlk_re = real_part(mat_Q[m][l][k]);
	else{
		Q_mlk_re = mat_C_re[m][l][k]/mat_B[m][l][k] + mat_d_re[l][k]*mat_beta[m][l][k]/mat_AC[l][k];
		Q_mlk_re *= -2;
	}
	return Q_mlk_re;
}
ex fn_Q_conj (int m, int l, int k){
	ex Q_mlk_conj;
	 /**  Conj(C_mlk) Conj(d_lk)
	 ** Conj(Q_mlk) = -2 ( ---------- + ------------- beta_mlk )
	 **  B_mlk AC_lk
	 **/
	if (my_is_zero(mat_B[m][l][k])==1 || my_is_zero(mat_AC[l][k])==1){
		stringstream err_msg;
		err_msg << "Q_" << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	if (is_a<numeric>(mat_Q[m][l][k]))
		Q_mlk_conj = mat_Q[m][l][k].conjugate();
	else{
		Q_mlk_conj = mat_C_conj[m][l][k]/mat_B[m][l][k] + mat_d_conj[l][k]*mat_beta[m][l][k]/mat_AC[l][k];
		Q_mlk_conj *= -2;
	}
	return Q_mlk_conj;
}

ex fn_P (int m, int l, int k){
	ex P_mlk;
	 /** A_mlk 
	 ** P_mlk = -2 [( ------- - alpha_lk) * (1 + beta_mlk.phi_mlk) - D_mlk . beta_mlk - phi_mlk ]
	 ** B_mlk
	 **/
	if (my_is_zero(mat_B[m][l][k])==1){
		stringstream err_msg;
		err_msg << "P_" << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	P_mlk = -(
	           (mat_A[m][l][k] / mat_B[m][l][k] - mat_alpha[l][k]) *(1 + mat_beta[m][l][k]*mat_phi[m][l][k])
	           - mat_D[m][l][k]*mat_beta[m][l][k] 
	           - mat_phi[m][l][k]
	         ) * 2;
	return P_mlk;
}

ex fn_E (int m, int l, int k){
	ex E_mlk;
	 /** d_lk C_mlk
	 ** E_mlk = -2 ( ------- + ------- phi_mlk )
	 ** AC_lk B_mlk
	 **/
	if (my_is_zero(mat_B[m][l][k])==1 || my_is_zero(mat_AC[l][k])==1){
		stringstream err_msg;
		err_msg << "E_" << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	E_mlk = -(
	           mat_d[l][k]/mat_AC[l][k] 
	           + mat_phi[m][l][k] * mat_C[m][l][k]/mat_B[m][l][k]
	         ) * 2;
	return E_mlk;
}
ex fn_E_im (int m, int l, int k){
	ex E_mlk_im;
	 /**  imag(d_lk) imag(C_mlk)
	 ** imag(E_mlk) = -2 ( ----------- + ------------- phi_mlk )
	 **  AC_lk B_mlk
	 **/
	if (my_is_zero(mat_B[m][l][k])==1 || my_is_zero(mat_AC[l][k])==1){
		stringstream err_msg;
		err_msg << "E_" << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	if (is_a<numeric>(mat_E[m][l][k]))
		E_mlk_im = imag_part(mat_E[m][l][k]);
	else{
		E_mlk_im = mat_d_im[l][k]/mat_AC[l][k] + mat_C_im[m][l][k] * mat_phi[m][l][k]/mat_B[m][l][k];
		E_mlk_im *= -2;
	}
	return E_mlk_im;
}
ex fn_E_re (int m, int l, int k){
	ex E_mlk_re;
	 /**  real(d_lk) real(C_mlk)
	 ** real(E_mlk) = -2 ( ----------- + ------------- phi_mlk )
	 **  AC_lk B_mlk
	 **/
	if (my_is_zero(mat_B[m][l][k])==1 || my_is_zero(mat_AC[l][k])==1){
		stringstream err_msg;
		err_msg << "E_" << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	if (is_a<numeric>(mat_E[m][l][k]))
		E_mlk_re = real_part(mat_E[m][l][k]);
	else{
		E_mlk_re = mat_d_re[l][k]/mat_AC[l][k] + mat_C_re[m][l][k] * mat_phi[m][l][k]/mat_B[m][l][k];
		E_mlk_re *= -2;
	}
	return E_mlk_re;
}

ex fn_F (int n, int m, int l, int k){
	ex F_nmlk;
	 /** C_nlk.B_mlk - C_mlk.B_nlk
	 ** F_nmlk = ---------------------------
	 ** A_nlk.B_mlk - A_mlk.B_nlk
	 **/
	ex denom = mat_A[n][l][k] * mat_B[m][l][k] - mat_A[m][l][k] * mat_B[n][l][k];
	if (my_is_zero(denom)==1){
		stringstream err_msg;
		err_msg << "F_" << n+1 << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	F_nmlk = mat_C[n][l][k]*mat_B[m][l][k] - mat_C[m][l][k]*mat_B[n][l][k];
	F_nmlk /= denom;
	return F_nmlk;
}
ex fn_F_im (int n, int m, int l, int k){
	ex F_nmlk_im;
	 /** Im(C_nlk).B_mlk - Im(C_mlk).B_nlk
	 ** Im(F_nmlk) = ------------------------------------
	 **  A_nlk.B_mlk - A_mlk.B_nlk
	 **/
	if (is_a<numeric>(mat_F[n][m][l][k]))
		F_nmlk_im = imag_part(mat_F[n][m][l][k]);
	else{
		ex denom = mat_A[n][l][k] * mat_B[m][l][k] - mat_A[m][l][k] * mat_B[n][l][k];
		if (my_is_zero(denom)==1){
			stringstream err_msg;
			err_msg << "F_" << n+1 << m+1 << l+1 << k+1 << " denominator = 0";
			throw std::runtime_error(err_msg.str());
		}
		F_nmlk_im = mat_C_im[n][l][k]*mat_B[m][l][k] - mat_C_im[m][l][k]*mat_B[n][l][k];
		F_nmlk_im /= denom;
	}
	return F_nmlk_im;
}
	//	II.6	Level 6 Variable Functions
ex fn_A0 (int m, int l, int k){
	 /**
	 ** A0_mlk = P_mlk * Im(E_mlk)
	 **/
	ex A0_mlk;
	A0_mlk = mat_P[m][l][k]*mat_E_im[m][l][k];
	return A0_mlk;
}

ex fn_B0 (int m, int l, int k){
	ex B0_mlk;
	 /**
	 ** B0_mlk = P_mlk.Im(m_k ^2) + rho.P_mlk + Re(Q_mlk).Im(E_mlk) - Im(Q_mlk).Re(E_mlk)
	 **/
	B0_mlk = mat_P[m][l][k]* imag_part(mat_msquare[k]) 
		+ Rho/*This is Feynman's Prescription*/ * mat_P[m][l][k] 
		+ mat_Q_re[m][l][k]*mat_E_im[m][l][k] 
		- mat_Q_im[m][l][k]*mat_E_re[m][l][k];
	return B0_mlk;
}

ex fn_C0 (int m, int l, int k){
	ex C0_mlk;
	 /**
	 ** C0_mlk = Im(Q_mlk) .Re(m_k ^2) + Re(Q_mlk).[Im(m_k ^2) + rho]
	 **/
	C0_mlk = mat_Q_im[m][l][k]*real_part(mat_msquare[k]) 
		+ mat_Q_re[m][l][k] * ( imag_part(mat_msquare[k]) + Rho /*This is Feynman's prescription*/);
	return C0_mlk;
}

ex fn_z1phi (int m, int l, int k){
	ex z1phi_mlk;
	 /** (E_mlk - Q_mlk.phi_mlk) + Sqrt( (E_mlk - Q_mlk.phi_mlk)^2 - 4.P_mlk.phi_mlk.(m_k ^2 - i.rho))
	 ** z1phi_mlk = ---------------------------------------------------------------------------------------------
	 **    2.P_mlk.phi_mlk
	 **/
	ex denom = mat_P[m][l][k]*mat_phi[m][l][k]*2;
	
	if (my_is_zero(denom)==1){
		stringstream err_msg;
		err_msg << "z1phi_" << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	
	z1phi_mlk = mat_E[m][l][k] - mat_Q[m][l][k]*mat_phi[m][l][k];
	z1phi_mlk += sqrt(
	                   pow(z1phi_mlk, 2)
	                   - mat_P[m][l][k]*mat_phi[m][l][k]*(mat_msquare[k] - I*Rho /*This is Feynman's prescription*/)*4 
	                 );
	z1phi_mlk /= mat_P[m][l][k] * mat_phi[m][l][k] * 2;
	
	return z1phi_mlk;
}
ex fn_z1phi_im (int m, int l, int k){
	ex z1phi_mlk_im;
	 /** (E_mlk - Q_mlk.phi_mlk) + Sqrt( (E_mlk - Q_mlk.phi_mlk)^2 - 4.P_mlk.phi_mlk.(m_k ^2 - i.rho))
	 ** z1phi_mlk = ---------------------------------------------------------------------------------------------
	 **    2.P_mlk.phi_mlk
	 **/
	ex denom = mat_P[m][l][k]*mat_phi[m][l][k] * 2;
	if (my_is_zero(denom)==1){
		stringstream err_msg;
		err_msg << "z1phi_" << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	z1phi_mlk_im = imag_part(mat_z1phi[m][l][k]);
	return z1phi_mlk_im;
}

ex fn_z2phi (int m, int l, int k){
	ex z2phi_mlk;
	 /** (E_mlk - Q_mlk.phi_mlk) - Sqrt( (E_mlk - Q_mlk.phi_mlk)^2 - 4.P_mlk.phi_mlk.(m_k ^2 - i.rho))
	 ** z2phi_mlk = ---------------------------------------------------------------------------------------------
	 **    2.P_mlk.phi_mlk
	 **/
	ex denom = 2*mat_P[m][l][k]*mat_phi[m][l][k];
	
	if (my_is_zero(denom)==1){
		stringstream err_msg;
		err_msg << "z2phi_" << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	
	z2phi_mlk = mat_E[m][l][k] - mat_Q[m][l][k]*mat_phi[m][l][k];
	z2phi_mlk -= sqrt(
	                   pow(z2phi_mlk,2)
	                   - mat_P[m][l][k]*mat_phi[m][l][k]*(mat_msquare[k] - I*Rho/*This is Feynman's prescrition*/)*4
	                 );
	z2phi_mlk /= mat_P[m][l][k]* mat_phi[m][l][k] * 2;
	
	return z2phi_mlk;
}
ex fn_z2phi_im (int m, int l, int k){
	ex z2phi_mlk_im;
	ex denom = mat_P[m][l][k]*mat_phi[m][l][k] * 2;
	
	if (my_is_zero(denom)==1){
		stringstream err_msg;
		err_msg << "z2phi_" << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	z2phi_mlk_im = imag_part(mat_z2phi[m][l][k]);
	return z2phi_mlk_im;
}

ex fn_z1beta (int m, int l, int k){
	ex z1beta_mlk;
	 /** (E_mlk - Q_mlk/beta_mlk) + Sqrt( (E_mlk - Q_mlk/beta_mlk)^2 - 4.P_mlk.phi_mlk.(m_k ^2 - i.rho)) ** z1beta_mlk = --------------------------------------------------------------------------------------------
	 **    2.P_mlk/beta_mlk
	 **/
	if (my_is_zero(mat_beta[m][l][k])==1 || my_is_zero(mat_P[m][l][k])==1){
		stringstream err_msg;
		err_msg << "z1beta_" << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	
	z1beta_mlk = mat_E[m][l][k] - mat_Q[m][l][k]/mat_beta[m][l][k];
	z1beta_mlk += sqrt(
	                    pow(z1beta_mlk, 2)
	                    - mat_P[m][l][k]*(mat_msquare[k] - I*Rho/*This is Feynman's prescription*/)/mat_beta[m][l][k] * 4
	                  );
	z1beta_mlk /= mat_P[m][l][k] / mat_beta[m][l][k] * 2;
	return z1beta_mlk;
}
ex fn_z1beta_im (int m, int l, int k){
	ex z1beta_mlk_im;
	ex denom = mat_P[m][l][k]/mat_beta[m][l][k] * 2;
	
	if (my_is_zero(mat_beta[m][l][k])==1 || my_is_zero(mat_P[m][l][k])==1){
		stringstream err_msg;
		err_msg << "z1beta_" << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	
	z1beta_mlk_im = imag_part(mat_z1beta[m][l][k]);
	return z1beta_mlk_im;
}

ex fn_z2beta (int m, int l, int k){
	ex z2beta_mlk;
	 /** (E_mlk - Q_mlk/beta_mlk) - Sqrt( (E_mlk - Q_mlk/beta_mlk)^2 - 4.P_mlk.phi_mlk.(m_k ^2 - i.rho)) ** z2beta_mlk = --------------------------------------------------------------------------------------------
	 **    2.P_mlk/beta_mlk
	 **/
	if (my_is_zero(mat_beta[m][l][k])==1 || my_is_zero(mat_P[m][l][k])==1){
		stringstream err_msg;
		err_msg << "z2beta_" << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	
	z2beta_mlk = mat_E[m][l][k] - mat_Q[m][l][k]/mat_beta[m][l][k];
	z2beta_mlk -= sqrt(
	                    pow(z2beta_mlk,2) 
	                    - mat_P[m][l][k]*(mat_msquare[k] - I*Rho/*This is Feynman's prescription*/)/mat_beta[m][l][k] * 4
	                  );
	z2beta_mlk /= mat_P[m][l][k] / mat_beta[m][l][k] * 2;
	
	return z2beta_mlk;
}
ex fn_z2beta_im (int m, int l, int k){
	ex z2beta_mlk_im;
	ex denom = mat_P[m][l][k]/mat_beta[m][l][k] * 2;
	
	if (my_is_zero(mat_beta[m][l][k])==1 || my_is_zero(mat_P[m][l][k])==1){
		stringstream err_msg;
		err_msg << "z2beta_" << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	z2beta_mlk_im = imag_part(mat_z2beta[m][l][k]);
	return z2beta_mlk_im;
}

ex fn_T1 (int n, int m, int l, int k){
	ex T1_nmlk;
	 /** Q_mlk + P_mlk.F_nmlk - beta_mlk.E_mlk)
	 ** T1_nmlk = ---------------------------------------- +
	 **  - 2.P_mlk
	 ** sqrt((Q_mlk+P_mlk.F_nmlk-beta_mlk.E_mlk)^2 - 4.P_mlk(Q_mlk.F_nmlk+beta_mlk.m_k^2-i.beta_mlk.rho))
	 ** + -------------------------------------------------------------------------------------------------
	 **    - 2.P_mlk
	 **/
	if (my_is_zero(mat_P[m][l][k])==1){
		stringstream err_msg;
		err_msg << "T1_" << n+1 << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	T1_nmlk = mat_Q[m][l][k] + mat_P[m][l][k]*mat_F[n][m][l][k] - mat_beta[m][l][k]*mat_E[m][l][k];
	T1_nmlk += sqrt(
	                 T1_nmlk*T1_nmlk 
	                 - mat_P[m][l][k]*(
	                                    mat_Q[m][l][k]*mat_F[n][m][l][k] 
	                                    + mat_beta[m][l][k]*mat_msquare[k] 
	                                    - I*mat_beta[m][l][k]*Rho /*This is Feynman's prescription*/
	                                  ) * 4
	               );
	T1_nmlk /= -mat_P[m][l][k] * 2;
	
	return T1_nmlk;
}

ex fn_T2 (int n, int m, int l, int k){
	ex T2_nmlk;
	 /** Q_mlk + P_mlk.F_nmlk - beta_mlk.E_mlk)
	 ** T2_nmlk = ---------------------------------------- +
	 **  - 2.P_mlk
	 ** sqrt((Q_mlk+P_mlk.F_nmlk-beta_mlk.E_mlk)^2 - 4.P_mlk(Q_mlk.F_nmlk+beta_mlk.m_k^2-i.beta_mlk.rho))
	 ** - -------------------------------------------------------------------------------------------------
	 **    - 2.P_mlk
	 **/
	if (my_is_zero(mat_P[m][l][k])==1){
		stringstream err_msg;
		err_msg << "T2_" << n+1 << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	T2_nmlk = mat_Q[m][l][k] + mat_P[m][l][k]*mat_F[n][m][l][k] - mat_beta[m][l][k]*mat_E[m][l][k];
	T2_nmlk -= sqrt(
	                 T2_nmlk*T2_nmlk 
	                 - 4*mat_P[m][l][k]*(
	                                      mat_Q[m][l][k]*mat_F[n][m][l][k] 
	                                      + mat_beta[m][l][k]*mat_msquare[k] 
	                                      - I*mat_beta[m][l][k]*Rho /*This is Feynman's prescription*/
	                                    ) 
	               );
	T2_nmlk /= -mat_P[m][l][k] * 2;
	return T2_nmlk;
}


	//	II.7	Level 7 Variable Functions
ex fn_OPlus (int n, int m, int l, int k){
	ex OPlus_nmlk;
	 /**   F_nmlk  1 - beta_mlk.phi_mlk
	 ** OPlus_nmlk = (f_lk.g_mlk+fminus_lk.g_mlk).ln(----------) - f_lk.g_mlk.ln(----------------------) +
	 **   beta_mlk  beta_mlk
	 **  1 - beta_mlk.phi_mlk
	 ** - f_lk.gminus_mlk.ln(----------------------) - 2pi.i.f_lk.g_mlk.Sign(Im(F_nmlk)) +
	 **  - beta_mlk
	 **   -P_mlk    Q_mlk
	 ** - (f_lk.g_mlk + fminus_lk.g_mlk)[ ln(----------) - ln(P_mlk) + 2.i.pi.theta(-P_mlk).Sign(Im(-------))
	 **   beta_mlk   P_mlk
	 **    P_mlk
	 **   + 2i.pi.theta(----------) + eta(z-z1beta_mlk,z-z2beta_mlk)] + 
	 **    beta_mlk
	 **     Q_mlk
	 ** + f_lk.g_mlk[ln(-P_mlk.phi_mlk) - ln(P_mlk) + 2i.pi.theta(-P_mlk).Sign(Im-------) +
	 **     P_mlk
	 **   + 2i.pi.theta(P_mlk.phi_mlk) + eta(z-z1phi_mlk, z-z2phi_mlk)] + 
	 **      Q_mlk
	 ** + f_lk.gminus_mlk[ln(-P_mlk.phi_mlk) - ln(-P_mlk) + 2i.pi.theta(P_mlk).Sign(Im-------) +
	 **      P_mlk
	 **   + 2i.pi.theta(P_mlk.phi_mlk) + eta(z-z1phi_mlk, z-z2phi_mlk)]
	 **/
	if (my_is_zero(mat_beta[m][l][k])==1 || my_is_zero(mat_P[m][l][k])==1){
		stringstream err_msg;
		err_msg << "OPlus_" << n+1 << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	ex term1, term2, term3, term4;
	
	term1 = (mat_f[l][k]*mat_g[m][l][k] +mat_fminus[l][k]*mat_g[m][l][k]) * log(mat_F[n][m][l][k]/mat_beta[m][l][k])
		- mat_f[l][k]*mat_g[m][l][k] * log( (-mat_beta[m][l][k]*mat_phi[m][l][k] + 1)/mat_beta[m][l][k]) 
		- mat_f[l][k]*mat_gminus[m][l][k] * log( (-mat_beta[m][l][k]*mat_phi[m][l][k] + 1)/(-mat_beta[m][l][k]) )
		+ Pi*I*mat_f[l][k]*mat_g[m][l][k] * my_csgn(mat_F_im[n][m][l][k]) * 2;
	
	term2 = -(mat_f[l][k]*mat_g[m][l][k] + mat_fminus[l][k]*mat_g[m][l][k]) * (
		log(-mat_P[m][l][k] / mat_beta[m][l][k]) 
		- log(mat_P[m][l][k]) 
		+ I*Pi*my_step(-mat_P[m][l][k])*my_csgn(imag_part(mat_Q[m][l][k] / mat_P[m][l][k])) * 2
		+ I*Pi*my_step(mat_P[m][l][k] / mat_beta[m][l][k]) * 2
		+ fn_eta_beta(m, l, k)
		);
	
	term3 = (mat_f[l][k]*mat_g[m][l][k])*(
	                                       log(-mat_P[m][l][k]*mat_phi[m][l][k]) 
	                                       - log(mat_P[m][l][k]) 
	                                       + my_step(-mat_P[m][l][k])*my_csgn(imag_part(mat_Q[m][l][k] / mat_P[m][l][k])) *2*I*Pi + my_step(mat_P[m][l][k]*mat_phi[m][l][k]) * 2*I*Pi
	                                       + fn_eta_phi(m,l,k)
	                                     );
	
	term4 = mat_f[l][k]*mat_gminus[m][l][k] * (
	                                            log(-mat_P[m][l][k]*mat_phi[m][l][k])
	                                            - log(-mat_P[m][l][k]) 
	                                            + my_step(mat_P[m][l][k])*my_csgn(imag_part(mat_Q[m][l][k] / mat_P[m][l][k])) *2*I*Pi
	                                            + my_step(mat_P[m][l][k]*mat_phi[m][l][k]) *2*I*Pi 
	                                            + fn_eta_phi(m, l, k)
	                                          );
	OPlus_nmlk = term1 + term2 + term3 + term4;	
	return OPlus_nmlk;
}

ex fn_OMinus (int n, int m, int l, int k){
	ex OMinus_nmlk;
	 /**   F_nmlk  -F_nmlk
	 ** OMinus_nmlk = -fminus_lk.gminus_mlk.ln(----------) - f_lk.gminus_mlk.ln(----------) +
	 **   beta_mlk  beta_mlk
	
	 ** - 2pi.i(fminus_lk.gminus_mlk+fminus_lk.g_mlk).Sign(Im(F_nmlk)) +
	 **    1 - beta_mlk.phi_mlk
	 ** + (fminus_lk.gminus_mlk + fminus_lk.g_mlk).ln(----------------------) +
	 **    beta_mlk
	 **   -P_mlk   Q_mlk
	 ** + fminus_lk.gminus_mlk.[ln(-----------) - ln(P_mlk) + 2i.pi.theta(-P_mlk).Sign(Im-------) +
	 **   beta_mlk   P_mlk
	 **    P_mlk
	 **   + 2i.pi.theta(----------) + eta(z-z1beta_mlk, z-z2beta_mlk) ] +
	 **    beta_mlk
	 **
	 **  -P_mlk   Q_mlk
	 ** + f_lk.gminus_mlk.[ln(-----------) - ln(-P_mlk) + 2i.pi.theta(P_mlk).Sign(Im-------) +
	 **  beta_mlk   P_mlk
	 **    P_mlk
	 **   + 2i.pi.theta(----------) + eta(z-z1beta_mlk, z-z2beta_mlk) ] +
	 **    beta_mlk
	 **
	 **       Q_mlk
	 ** + (fminus_lk.gminus_mlk + fminus_lk.g_mlk).[ln(-P_mlk.phi_mlk) - ln(P_mlk) + 2i.pi.theta(-P_mlk).Sign(Im-------) +
	 **       P_mlk
	 **
	 **   + 2i.pi.theta(P_mlk.phi_mlk) + eta(z-z1phi_mlk, z-z2phi_mlk) ]
	 **
	 **/
	if (my_is_zero(mat_beta[m][l][k])==1 || my_is_zero(mat_P[m][l][k])==1){
		stringstream err_msg;
		err_msg << "OMinus_" << n+1 << m+1 << l+1 << k+1 << " denominator = 0";
		throw std::runtime_error(err_msg.str());
	}
	ex term1, term2, term3, term4;
	
	term1 = -mat_fminus[l][k]*mat_gminus[m][l][k]*log(mat_F[n][m][l][k] / mat_beta[m][l][k])
		- mat_f[l][k]*mat_gminus[m][l][k]*log(-mat_F[n][m][l][k]/mat_beta[m][l][k])
		- (mat_fminus[l][k]*mat_gminus[m][l][k] + mat_fminus[l][k]*mat_g[m][l][k])*my_csgn(mat_F_im[n][m][l][k]) *2*I*Pi 
		+ (mat_fminus[l][k]*mat_gminus[m][l][k] + mat_fminus[l][k]*mat_g[m][l][k])*log((1-mat_beta[m][l][k]*mat_phi[m][l][k])/mat_beta[m][l][k]);
	
	term2 = mat_fminus[l][k]*mat_gminus[m][l][k] * (
	                                                 log(-mat_P[m][l][k]/mat_beta[m][l][k]) 
	                                                 - log(mat_P[m][l][k])
	                                                 + my_step(-mat_P[m][l][k])*my_csgn(imag_part(mat_Q[m][l][k]/mat_P[m][l][k])) *2*I*Pi 
	                                                 + my_step(mat_P[m][l][k]/mat_beta[m][l][k]) *2*I*Pi 
	                                                 + fn_eta_beta(m, l, k)
	                                               );
	
	term3 = mat_f[l][k]*mat_gminus[m][l][k]*(
	                                          log(-mat_P[m][l][k]/mat_beta[m][l][k])
	                                          - log(-mat_P[m][l][k])
	                                          + my_step(mat_P[m][l][k])*my_csgn(imag_part(mat_Q[m][l][k]/mat_P[m][l][k])) *2*I*Pi 
	                                          + my_step(mat_P[m][l][k]/mat_beta[m][l][k]) *2*I*Pi 
	                                          + fn_eta_beta(m, l, k)
	                                        );
	
	term4 = -(mat_fminus[l][k]*mat_gminus[m][l][k] + mat_fminus[l][k]*mat_g[m][l][k])*(
		log(-mat_P[m][l][k]*mat_phi[m][l][k])
		- log(mat_P[m][l][k]) 
		+ my_step(-mat_P[m][l][k])*my_csgn(imag_part(mat_Q[m][l][k]/mat_P[m][l][k])) *2*I*Pi 
		+ my_step(mat_P[m][l][k]*mat_phi[m][l][k]) *2*I*Pi 
		+ fn_eta_phi(m, l, k)
		);
	OMinus_nmlk = term1 + term2 + term3 + term4;
	return OMinus_nmlk;
}


	//	II.8	Level 8 Function
void fn_1Loop4Pt_calc(){
	int k, l, m, n;
	for(k=0; k<4; k++){
		for(l=0; l<4; l++){
			mat_a[k][l] = 0;
			mat_b[k][l] = 0;
			mat_c[k][l] = 0;
			mat_d[k][l] = 0; mat_d_im[k][l] = 0; mat_d_re[k][l] = 0; mat_d_conj[k][l] = 0;
			mat_AC[k][l] = 0;
			mat_alpha[k][l] = 0;
			mat_f[k][l] = 0;
			mat_fminus[k][l] = 0;
			for(m=0; m<4; m++)
			{
				mat_A[k][l][m] = 0;
				mat_B[k][l][m] = 0;
				mat_C[k][l][m] = 0; mat_C_im[k][l][m] = 0; mat_C_re[k][l][m] = 0; mat_C_conj[k][l][m] = 0;
				mat_D[k][l][m] = 0;
				mat_beta[k][l][m] = 0;
				mat_phi[k][l][m] = 0;
				mat_g[k][l][m] = 0;
				mat_gminus[k][l][m] = 0;
				mat_Q[k][l][m] = 0; mat_Q_im[k][l][m] = 0; mat_Q_re[k][l][m] = 0; mat_Q_conj[k][l][m] = 0;	
				mat_P[k][l][m] = 0;
				mat_E[k][l][m] = 0; mat_E_im[k][l][m] = 0; mat_E_re[k][l][m] = 0;
				
				mat_A0[k][l][m] = 0;
				mat_B0[k][l][m] = 0;
				mat_C0[k][l][m] = 0;
				mat_z1phi[k][l][m] = 0;mat_z1phi_im[k][l][m] = 0;
				mat_z2phi[k][l][m] = 0; mat_z2phi_im[k][l][m] = 0;
				mat_z1beta[k][l][m] = 0; mat_z1beta_im[k][l][m] = 0;
				mat_z2beta[k][l][m] = 0; mat_z2beta_im[k][l][m] = 0; 
				for(n=0;n<4;n++){
					mat_F[k][l][m][n] = 0; mat_F_im[k][l][m][n] = 0;
					mat_T1[k][l][m][n] = 0;
					mat_T2[k][l][m][n] = 0;
					mat_OPlus[k][l][m][n] = 0;
					mat_OMinus[k][l][m][n] = 0;
				}
			}
		}
	}
	
	///*Calculate Level 1 Terms*/
	for(k = 0; k<4; k++) for(l = 0; l<4; l++) if (l!=k){
		mat_a[l][k] = fn_a(l, k);
		mat_b[l][k] = fn_b(l, k);
		mat_c[l][k] = fn_c(l, k);
		mat_d[l][k] = fn_d(l, k); mat_d_im[l][k] = fn_d_im(l, k); mat_d_re[l][k] = fn_d_re(l, k); mat_d_conj[l][k] = fn_d_conj(l, k);
	}
	///*Calculate Level 2 Terms*/
	for(k = 0; k<4; k++) for(l = 0; l<4; l++) if (l!=k){
		mat_AC[l][k] = fn_AC(l, k);
		mat_alpha[l][k] = fn_alpha(l, k);
	}
	
	///*Calculate Level 3 Terms*/
	for(k = 0; k<4; k++) for(l = 0; l<4; l++) if (l!=k){
		mat_f[l][k] = fn_f(l, k);
		mat_fminus[l][k] = fn_fminus(l,k);
	}
	for(k = 0; k<4; k++) for(l = 0; l<4; l++) if (l!=k) for(m = 0; m<4; m++) if(m!=l && m!=k){
		mat_A[m][l][k] = fn_A(m, l, k);
		mat_B[m][l][k] = fn_B(m, l, k);
		mat_C[m][l][k] = fn_C(m, l, k); mat_C_im[m][l][k] = fn_C_im(m, l, k); mat_C_re[m][l][k] = fn_C_re(m, l, k); mat_C_conj[m][l][k] = fn_C_conj(m, l, k);
		mat_D[m][l][k] = fn_D(m, l, k);
	}
	
	///*Calculate Level 4 Terms*/
	for(k = 0; k<4; k++) for(l = 0; l<4; l++) if (l!=k) for(m = 0; m<4; m++) if(m!=l && m!=k){
		mat_beta[m][l][k]	=fn_beta(m, l, k);
		mat_phi[m][l][k]	=fn_phi(m, l, k);
		mat_g[m][l][k]		=fn_g(m, l, k);
		mat_gminus[m][l][k]	=fn_gminus(m, l, k);
	}
	///*Calculate Level 5 Terms*/
	for(k=0; k<4; k++) for(l=0; l<4; l++) if (l!=k) for(m = 0; m<4; m++) if(m!=l && m!=k){
		mat_Q[m][l][k]	=fn_Q(m, l, k); mat_Q_im[m][l][k]=fn_Q_im(m, l, k); mat_Q_re[m][l][k]=fn_Q_re(m, l, k); mat_Q_conj[m][l][k]=fn_Q_conj(m, l, k);
		mat_P[m][l][k]	=fn_P(m, l, k);
		mat_E[m][l][k]	=fn_E(m, l, k); mat_E_im[m][l][k]=fn_E_im(m, l, k); mat_E_re[m][l][k]=fn_E_re(m, l, k);
	}
	for(k=0; k<4; k++) for(l=0; l<4; l++) if (l!=k) for(m=0; m<4; m++) if(m!=l && m!=k) for(n=0; n<4; n++) if(n!=m && n!=l && n!=k){
		mat_F[n][m][l][k] =fn_F(n, m, l, k); mat_F_im[n][m][l][k] =fn_F_im(n, m, l, k);
	}
	// We skip it here, for calculating using Vegas
//blah blah blah
	D0 = 0;
	return;
//blah blah blah
		///*Calculate Level 6 Terms*/
	for(k = 0; k<4; k++) for(l = 0; l<4; l++) if (l!=k) for(m = 0; m<4; m++) if(m!=l && m!=k){
		mat_A0[m][l][k]	=fn_A0(m, l, k);
		mat_B0[m][l][k]	=fn_B0(m, l, k);
		mat_C0[m][l][k]	=fn_C0(m, l, k);
		mat_z1phi[m][l][k]	=fn_z1phi(m, l, k); mat_z1phi_im[m][l][k] =fn_z1phi_im(m, l, k);
		mat_z2phi[m][l][k]	=fn_z2phi(m, l, k); mat_z2phi_im[m][l][k] =fn_z2phi_im(m, l, k);
		mat_z1beta[m][l][k]	=fn_z1beta(m, l, k); mat_z1beta_im[m][l][k] =fn_z1beta_im(m, l, k);
		mat_z2beta[m][l][k]	=fn_z2beta(m, l, k); mat_z2beta_im[m][l][k] =fn_z2beta_im(m, l, k);
	}
	for(k = 0; k<4; k++) for(l = 0; l<4; l++) if (l!=k) for(m = 0; m<4; m++) if(m!=l && m!=k) for(n = 0; n<4; n++) if(n!=m && n!=l && n!=k){
		mat_T1[n][m][l][k] = fn_T1(n, m, l, k);
		mat_T2[n][m][l][k] = fn_T2(n, m, l, k);
	}
	
	///*Calculate Level 7 Terms*/
	for(k = 0; k<4; k++) for(l = 0; l<4; l++) if (l!=k) for(m = 0; m<4; m++) if(m!=l && m!=k) for(n = 0; n<4; n++) if(n!=m && n!=l && n!=k){
		mat_OPlus[n][m][l][k] = fn_OPlus(n, m, l, k);
		mat_OMinus[n][m][l][k] = fn_OMinus(n, m, l, k);
	}
	
	D0 = 0;
	ex term1 = 0, term2 = 0, 
		term2_1_2 = 0, term2_3_4 = 0, term2_5_7 = 0, term2_5_7_1 = 0, term2_5_7_2 = 0, term2_6_8 = 0, term2_11_13 = 0, term2_12_14 = 0, term2_9_15 = 0, term2_10_16 = 0;
	for(k = 0; k<4; k++) for(l = 0; l<4; l++) if (l!=k) for(m = 0; m<4; m++) if(m!=l && m!=k) for(n = 0; n<4; n++) if(n!=m && n!=l && n!=k){
		term1 = mat_B[m][l][k]/(mat_B[m][l][k]*mat_A[n][l][k] - mat_B[n][l][k]*mat_A[m][l][k]);
		term1 *=(1-myfn_delta(mat_AC[l][k]))*(1-myfn_delta(mat_B[m][l][k]))*abs(1-mat_beta[m][l][k]*mat_phi[m][l][k])*(-1/mat_P[m][l][k]);
		
		term2_1_2 += /*1*/mat_OPlus[n][m][l][k]*fn_GPositive(mat_T1[n][m][l][k],mat_T2[n][m][l][k]);
		term2_1_2 += /*2*/mat_OMinus[n][m][l][k]*fn_GNegative(mat_T1[n][m][l][k],mat_T2[n][m][l][k]);
		term2_3_4 += /*3*/ mat_fminus[l][k]*mat_g[m][l][k]*my_step(mat_Q_im[m][l][k].conjugate())*fn_ThetaG(-mat_A0[m][l][k], -mat_B0[m][l][k], -mat_C0[m][l][k], mat_T1[n][m][l][k], mat_T2[n][m][l][k]) *2*I*Pi ;
		term2_3_4 +=/*4*/ -mat_f[l][k]*mat_gminus[m][l][k]*my_step(-mat_Q_im[m][l][k].conjugate())*fn_ThetaG(mat_A0[m][l][k], mat_B0[m][l][k], mat_C0[m][l][k], mat_T1[n][m][l][k], mat_T2[n][m][l][k]) *2*I*Pi ;
		term2_5_7 +=/*5*/ -(mat_f[l][k]*mat_g[m][l][k]+mat_fminus[l][k]*mat_g[m][l][k])*fn_LnGPositive(mat_T1[n][m][l][k], mat_T2[n][m][l][k], mat_z1beta[m][l][k]);
		term2_5_7_1 = -(mat_f[l][k]*mat_g[m][l][k]+mat_fminus[l][k]*mat_g[m][l][k])*fn_LnGPositive(mat_T1[n][m][l][k], mat_T2[n][m][l][k], mat_z1beta[m][l][k]);
		term2_6_8 +=/*6*/ -(mat_f[l][k]*mat_g[m][l][k]+mat_fminus[l][k]*mat_g[m][l][k])*fn_LnGPositive(mat_T1[n][m][l][k], mat_T2[n][m][l][k], mat_z2beta[m][l][k]);
		term2_5_7 +=/*7*/ (mat_f[l][k]*mat_g[m][l][k]+mat_f[l][k]*mat_gminus[m][l][k])*fn_LnGPositive(mat_T1[n][m][l][k], mat_T2[n][m][l][k], mat_z1phi[m][l][k]);
		term2_5_7_2 = (mat_f[l][k]*mat_g[m][l][k]+mat_f[l][k]*mat_gminus[m][l][k])*fn_LnGPositive(mat_T1[n][m][l][k], mat_T2[n][m][l][k], mat_z1phi[m][l][k]);
		term2_6_8 +=/*8*/ (mat_f[l][k]*mat_g[m][l][k]+mat_f[l][k]*mat_gminus[m][l][k])*fn_LnGPositive(mat_T1[n][m][l][k], mat_T2[n][m][l][k], mat_z2phi[m][l][k]);
		term2_9_15 +=/*9*/ (mat_fminus[l][k]*mat_g[m][l][k]-mat_f[l][k]*mat_gminus[m][l][k])*fn_LnGPositive(mat_T1[n][m][l][k], mat_T2[n][m][l][k], -mat_Q[m][l][k]/mat_P[m][l][k]);
		term2_10_16 +=/*10*/ -(mat_f[l][k]*mat_g[m][l][k]+mat_f[l][k]*mat_gminus[m][l][k])*fn_LnGPositive(mat_T1[n][m][l][k], mat_T2[n][m][l][k], -mat_F[n][m][l][k]/(1-mat_beta[m][l][k]*mat_phi[m][l][k]));
		term2_11_13 +=/*11*/ (mat_fminus[l][k]*mat_gminus[m][l][k]+mat_f[l][k]*mat_gminus[m][l][k])*fn_LnGNegative(mat_T1[n][m][l][k], mat_T2[n][m][l][k], mat_z1beta[m][l][k]);
		term2_12_14 +=/*12*/ (mat_fminus[l][k]*mat_gminus[m][l][k]+mat_f[l][k]*mat_gminus[m][l][k])*fn_LnGNegative(mat_T1[n][m][l][k], mat_T2[n][m][l][k], mat_z2beta[m][l][k]);
		term2_11_13 +=/*13*/ -(mat_fminus[l][k]*mat_gminus[m][l][k]+mat_fminus[l][k]*mat_g[m][l][k])*fn_LnGNegative(mat_T1[n][m][l][k], mat_T2[n][m][l][k], mat_z1phi[m][l][k]);
		term2_12_14 +=/*14*/ -(mat_fminus[l][k]*mat_gminus[m][l][k]+mat_fminus[l][k]*mat_g[m][l][k])*fn_LnGNegative(mat_T1[n][m][l][k], mat_T2[n][m][l][k], mat_z2phi[m][l][k]);
		term2_9_15 +=/*15*/ (mat_fminus[l][k]*mat_g[m][l][k]-mat_f[l][k]*mat_gminus[m][l][k])*fn_LnGNegative(mat_T1[n][m][l][k], mat_T2[n][m][l][k], -mat_Q[m][l][k]/mat_P[m][l][k]);
		term2_10_16 +=/*16*/ (mat_fminus[l][k]*mat_gminus[m][l][k]+mat_fminus[l][k]*mat_g[m][l][k])*fn_LnGNegative(mat_T1[n][m][l][k], mat_T2[n][m][l][k], -mat_F[n][m][l][k]/(1-mat_beta[m][l][k]*mat_phi[m][l][k]));
		term2 = term2_1_2 + term2_3_4 + term2_5_7 + term2_6_8 + term2_11_13 + term2_12_14 + term2_9_15 + term2_10_16;
		D0 += term1 * term2;
	}
	
	D0 *= I*pow(Pi,2);
}

	//	II.9	Misc. Functions
ex fn_delta_eval (const ex & x){	// Delta Chroenacker
	if (is_a<numeric>(x))
		return my_is_zero(x);
	else
		return myfn_delta(x).hold();
}

REGISTER_FUNCTION(myfn_delta, eval_func(fn_delta_eval));

ex myfn_IvZero_eval(const ex & x){
	if (my_is_zero(x) == 1)
		return 0;
	else if (my_is_zero(x) == 0)
		return 1;
	return myfn_IvZero(x).hold();
}

REGISTER_FUNCTION(myfn_IvZero, eval_func(myfn_IvZero_eval));

ex myfn_Zero_eval(const ex & x){
	if (is_a<numeric>(x))
		return my_is_zero(x);
	else
		return myfn_Zero(x).hold();
}

REGISTER_FUNCTION(myfn_Zero, eval_func(myfn_Zero_eval));

ex myfn_IvTheta_eval(const ex & x){
	if (is_a<numeric>(x)){
		if (real_part(x).info(info_flags::negative))
			return 1;
		else
			return 0;
	}
	else
		return myfn_IvTheta(x).hold();
}

REGISTER_FUNCTION(myfn_IvTheta, eval_func(myfn_IvTheta_eval));

ex myfn_Sign0_eval(const ex & x){
	if (is_a<numeric>(x)){
		if (real_part(x).info(info_flags::negative))
			return -1;
		else
			return 1;
	}
	else
		return myfn_Sign0(x).hold();
}

REGISTER_FUNCTION(myfn_Sign0, eval_func(myfn_Sign0_eval));


ex fn_eta_eval(const ex &x, const ex &y){
	 /** *****************************************************************************
	 **	eta(x,y) = 2pi.i.{
	 **			theta(-Im(x)).theta(-Im(y)).theta(Im(x.y))
	 **			- 
	 **			theta(Im(x)).theta(Im(y)).theta(-Im(x.y))
	 **			}
	 *********************************************************************************/
	if (is_a<numeric>(x) && is_a<numeric>(y)){
		ex eta;
		eta = 2*Pi*I*(
		               my_step(-imag_part(x))*my_step(-imag_part(y))*my_step(imag_part(x*y))
		               -
		               my_step(imag_part(x))*my_step(imag_part(y))*my_step(-imag_part(x*y))
		             );
		return eta;
	}
	return fn_eta(x,y).hold();
}
REGISTER_FUNCTION(fn_eta, eval_func(fn_eta_eval));

ex fn_eta_beta(int m, int l, int k){
	 /** *****************************************************************************
	 **	eta(z-z1beta_mlk,z-z2beta_mlk) = {
	 **			theta[Im(z1beta_mlk)].theta[Im(z2beta_mlk)].theta[-P_mlk/beta_mlk]
	 **			- 
	 **			theta[-Im(z1beta_mlk)].theta[-Im(z2beta_mlk)].theta[P_mlk/beta_mlk]
	 **			}
	 *********************************************************************************/
	ex eta;
	if(my_is_zero(mat_z1beta_im[m][l][k])==1)
		return 0;
	
	eta = my_step(mat_z1beta_im[m][l][k])*my_step(mat_z2beta_im[m][l][k])*my_step(-mat_P[m][l][k]/mat_beta[m][l][k])
		-
		my_step(-mat_z1beta_im[m][l][k])*my_step(-mat_z2beta_im[m][l][k])*my_step(mat_P[m][l][k]/mat_beta[m][l][k]);
	return eta *2*I*Pi ;
}

ex fn_eta_phi(int m, int l, int k){
	 /** *****************************************************************************
	 **	eta(z-z1phimlk,z-z2phi_mlk) = {
	 **			theta[Im(z1phi_mlk)].theta[Im(z2phi_mlk)].theta[-P_mlk*phi_mlk]
	 **			- 
	 **			theta[-Im(z1phi_mlk)].theta[-Im(z2phi_mlk)].theta[P_mlk*phi_mlk]
	 **			}
	 *********************************************************************************/
	ex eta;
	if(my_is_zero(mat_z1phi_im[m][l][k])==1)
		return 0;
	
	eta = my_step(mat_z1phi_im[m][l][k])*my_step(mat_z2phi_im[m][l][k])*my_step(-mat_P[m][l][k]*mat_phi[m][l][k])
		-
		my_step(-mat_z1phi_im[m][l][k])*my_step(-mat_z2phi_im[m][l][k])*my_step(mat_P[m][l][k]*mat_phi[m][l][k]);
	return eta *2*I*Pi ;
}

ex fn_R(const ex &x, const ex &y){
	 /** *****************************************************************************
	 **		ln(x) - ln(y)
	 **	R = ----------------------
	 **		 x - y
	 *********************************************************************************/
	if(x == y)
		throw std::runtime_error("R Function's Denominator = 0");
	
	return (log(x) - log(y)) / (x - y);
}

ex fn_GPositive(const ex &T1, const ex &T2){
	 /** *****************************************************************************
	 **	G+ = R(-T1, -T2)
	 *********************************************************************************/
	return fn_R(-T1, -T2);
}

ex fn_GNegative(const ex &T1, const ex &T2){
	 /** *****************************************************************************
	 **	G- = R(T1, T2)
	 *********************************************************************************/
	return fn_R(T1, T2);
}

ex fn_LnGPositive(const ex &T1, const ex &T2, const ex &z0){
	ex LnGPositive = 0;
	if(my_is_zero(T1)==1 || my_is_zero(T2)==1 || T1==T2)
		throw std::runtime_error("LnGPositive Function's Denominator = 0");
	ex log_mT1 = log(-T1), log_mT2 = log(-T2);
	
	LnGPositive=fn_R(-T1, -T2);
	LnGPositive += (
	                 log_mT1 - log_mT2 - pow(log_mT1,2)/2 + pow(log_mT2,2)/2 + Li2(-z0/T2 + 1) - Li2(-z0/T1 + 1) 
	               )/(T1-T2);
	LnGPositive += (
	                 log_mT2*( 
	                           fn_eta(T2-z0, (ex)1/(T2 + 1)) - fn_eta(T2-z0, (ex)1/T2)
	                         )
	                 -log_mT1*( 
	                            fn_eta(T1-z0, (ex)1/(T1+1)) - fn_eta(T1-z0, (ex)1/T1)
	                          )
	               )/(T1-T2);
	LnGPositive +=(
	                log(-z0/T2 + 1)*fn_eta(-z0, -(ex)1/T2) - log(-z0/T1 + 1)*fn_eta(-z0, -(ex)1/T1)
	              )/(T1-T2);
	return LnGPositive;
}

ex fn_LnGNegative(const ex &T1, const ex &T2, const ex &z0){
	ex LnGNegative;
	if(my_is_zero(T1)==1 || my_is_zero(T2)==1 || T1==T2)
		throw std::runtime_error("LnGPositive Function's Denominator = 0");
	
	ex log_T1 = log(T1), log_T2 = log(T2);
	LnGNegative=(-log(z0) + log(-z0) + 1)*fn_R(T1, T2);
	LnGNegative += - (
	                   log_T1 - log_T2 - (pow(log_T1,2))/2 + (pow(log_T2,2))/2 + Li2(-z0/T2 + 1) - Li2(-z0/T1 + 1)
	                 )/(T1-T2);
	LnGNegative += - (
	                   log_T2*( 
	                            fn_eta(z0-T2, (ex)1/(-T2 + 1)) - fn_eta(z0-T2, -(ex)1/T2)
	                          )
	                   -log_T1*(
	                             fn_eta(-T1+z0, (ex)1/(-T1 + 1)) - fn_eta(-T1+z0, -(ex)1/T1)
	                           )
	                 )/(T1-T2);
	LnGNegative += - (
	                   log(-z0/T2 + 1)*fn_eta(z0, (ex)1/T2) - log(-z0/T1 + 1)*fn_eta(z0, (ex)1/T1)
	                 )/(T1-T2);
	return LnGNegative;
}

ex fn_ThetaG(const ex &A0, const ex &B0, const ex &C0, const ex &T1, const ex &T2){
	ex ThetaG;
	ex m_delta, z1, z2; /*Note: delta is already a function, use m_delta instead*/
	// if(my_is_zero(B0))
	//	throw std::runtime_error("ThetaG Function's Denominator = 0");
	if(is_a<numeric>(A0) && is_a<numeric>(B0) && is_a<numeric>(C0) && is_a<numeric>(T1) && is_a<numeric>(T2)){
		//the input is numeric, return numeric
		if(my_is_zero(A0)==1 && my_is_zero(B0)==1){
			if(real_part(C0)>0 or my_is_zero(C0)==1) // Case 1
				return fn_R(T1, T2) + fn_R(-T1,-T2);
			else // case 2
				return 0;
		}
		else if(my_is_zero(A0)==1){
			if(real_part(B0).info(info_flags::positive)) // Case 3
				return fn_R(-C0/B0-T1,-C0/B0-T2);
			else // Case 4
				return fn_R(C0/B0 + T1, C0/B0 + T2);
		}
		m_delta = B0*B0 - A0*C0*4;
		z1 = (-B0 + sqrt(m_delta))/(A0*2);
		z2 = (-B0 - sqrt(m_delta))/(A0*2);
		if(real_part(A0).info(info_flags::positive)){
			if ( !real_part(m_delta).info(info_flags::positive)) // case 5
				return fn_R(T1, T2) + fn_R(-T1, -T2);
			else if (real_part(z2)>real_part(z1)) // case 6
				return fn_R(T1-z1, T2-z1) + fn_R(z2-T1, z2-T2);
			else // case 7
				return fn_R(T1-z2, T2-z2) + fn_R(z1-T1, z1-T2);
		}
		else{
			if ( !real_part(m_delta).info(info_flags::positive)) // case 8
				return 0;
			else if (real_part(z2)>real_part(z1)) // case 9
				return fn_R(z1-T1, z1-T2) - fn_R(z2-T1, z2-T2);
			else // case 10
				return fn_R(z2-T1, z2-T2) - fn_R(z1-T1, z1-T2);
		}
	}
	else{// the input is not numeric, return symbol
		m_delta = B0*B0 - A0*C0*4;
		
		z1 = ((-B0 + sqrt(m_delta))/(A0*2 + Rho2/*This is Khiem's prescription*/)) * myfn_IvZero(A0) * myfn_Zero(-m_delta) - (C0/(B0 + Rho2/*This is Khiem's prescription*/)) * myfn_Zero(A0) * myfn_IvZero(B0);
		z2 = ((-B0 - sqrt(m_delta))/(A0*2 + Rho2/*This is Khiem's prescription*/)) * myfn_IvZero(A0) * myfn_Zero(-m_delta) - (C0/(B0 + Rho2/*This is Khiem's prescription*/)) * myfn_Zero(A0) * myfn_IvZero(B0);
		
		ThetaG = - myfn_Zero(A0) * myfn_Zero(B0) * myfn_Zero(C0) + 1;
		ThetaG *= - myfn_Zero(A0) * my_step(-m_delta) + 1;
		ThetaG *= fn_R( (T1-z1)*myfn_Sign0(-A0), (T2-z1)*myfn_Sign0(-A0) )
			+ myfn_Sign0(A0) * fn_R((T1-z1)*myfn_Sign0(A0), (T2-z1)*myfn_Sign0(A0))
			- myfn_IvZero(B0) * myfn_Zero(A0) * fn_R( (T1 + C0/B0)*myfn_Sign0(B0), (T2 + C0/B0)* myfn_Sign0(B0));
		return ThetaG;
	}
}

//}// Namespace xloops

void input_manual(lst &q, lst &m, ex &Rho, ex &Rho1, ex &Rho2){
	float m1, m2, m3, m4, q10, q20, q21, q30, q31, q32;
	cout << "m1^2 = ";cin >> m1;
	cout << "m2^2 = ";cin >> m2;
	cout << "m3^2 = ";cin >> m3;
	cout << "m4^2 = ";cin >> m4;
	cout << "q10 = ";cin >> q10;
	cout << "q20 = ";cin >> q20;
	cout << "q21 = ";cin >> q21;
	cout << "q30 = ";cin >> q30;
	cout << "q31 = ";cin >> q31;
	cout << "q32 = ";cin >> q32;
	cout << "Rho = ";cin >> Rho;
	cout << "Rho1 = ";cin >> Rho1;
	cout << "Rho2 = ";cin >> Rho2;
	m = lst(m1, m2, m3, m4);
	q = lst(q10, q20, q21, q30, q31, q32);
}

void input_option(int argc, char *argv[], lst &q, lst &m, ex &Rho, ex &Rho1, ex &Rho2){
	/*Use optget_long to get the long options with their values*/
	ex m1, m2, m3, m4, q10, q20, q21, q30, q31, q32;
	
	m1 = symbol("m1"); m2 = symbol("m2"); m3 = symbol("m3"); m4 = symbol("m4");
	q10 = symbol("q10");
	q20 = symbol("q20"); q21 = symbol("q21");
	q30 = symbol("q30"); q31 = symbol("q31"); q32 = symbol("q32");
	Rho = 1e-25; Rho1 = 1e-24; Rho2 = 1e-23;
	
	// Store the P's for loop tool
	ex p1s=symbol("p1s"),p2s=symbol("p2s"),p3s=symbol("p3s"),p4s=symbol("p4s"),p12s=symbol("p12s"),p23s=symbol("p23s");
	
	
	mat_msquare[0]=symbol("m1");	mat_msquare[1]=symbol("m2");	mat_msquare[2]=symbol("m3");	mat_msquare[3]=symbol("m4");
	mat_q[0][0]=symbol("q10");mat_q[0][1]=0;mat_q[0][2]=0;mat_q[0][3]=0;
	mat_q[1][0]=symbol("q20");mat_q[1][1]=symbol("q21");mat_q[1][2]=0;mat_q[1][3]=0;
	mat_q[2][0]=symbol("q30");mat_q[2][1]=symbol("q31");mat_q[2][2]=symbol("q32");mat_q[2][3]=0;
	
	
	
	int c;
	int digit_optind = 0;
	double temp_numcast = 0;
	short int loopToolInput_Indicator = 0; // Whether input in P (as for loop tool):
	short int xloopsInput_Indicator = 0; // Whether input in Q (xloop type)
	
	
	while (1)
	{
		int this_option_optind = optind ? optind : 1;
		int option_index = 0;
		static struct option long_options[] =
		{
			{"m1s", 1, 0, 0},
			{"m2s", 1, 0, 0},
			{"m3s", 1, 0, 0},
			{"m4s", 1, 0, 0},
			{"q10", 1, 0, 0},
			{"q20", 1, 0, 0},
			{"q21", 1, 0, 0},
			{"q30", 1, 0, 0},
			{"q31", 1, 0, 0},
			{"q32", 1, 0, 0},
			{"p1s", 1, 0, 0},//P1 Square
			{"p2s", 1, 0, 0},//P2 Square
			{"p3s", 1, 0, 0},//P3 Square
			{"p4s", 1, 0, 0},//P4 Square
			{"p12s", 1, 0, 0},//(P1+P2) Square
			{"p23s", 1, 0, 0},//(P2+P3) Square
			{"rho", 1, 0, 0},
			{"rho1", 1, 0, 0},
			{"rho2", 1, 0, 0},
			{0, 0, 0, 0}
		};
		c = getopt_long (argc, argv, "",
		                 long_options, &option_index);
		if (c == -1)
			break;
		
		switch (c)
		{
		case 0:
			temp_numcast = strtod(optarg, NULL);
			switch(option_index)
			{
			case 0: m1 = temp_numcast;
				break;
			case 1: m2 = temp_numcast;
				break;
			case 2: m3 = temp_numcast;
				break;
			case 3: m4 = temp_numcast;
				break;
			case 4: q10 = temp_numcast;
				xloopsInput_Indicator = 1;
				break;
			case 5: q20 = temp_numcast;
				xloopsInput_Indicator = 1;
				break;
			case 6: q21 = temp_numcast;
				xloopsInput_Indicator = 1;
				break;
			case 7: q30 = temp_numcast;
				xloopsInput_Indicator = 1;
				break;
			case 8: q31 = temp_numcast;
				xloopsInput_Indicator = 1;
				break;
			case 9: q32 = temp_numcast;
				xloopsInput_Indicator = 1;
				break;
			case 10: p1s = temp_numcast; 
				loopToolInput_Indicator = 1;
				break;
			case 11: p2s = temp_numcast; 
				loopToolInput_Indicator = 1;
				break;
			case 12: p3s = temp_numcast; 
				loopToolInput_Indicator = 1;
				break;
			case 13: p4s = temp_numcast; 
				loopToolInput_Indicator = 1;
				break;
			case 14: p12s = temp_numcast; 
				loopToolInput_Indicator = 1;
				break;
			case 15: p23s = temp_numcast; 
				loopToolInput_Indicator = 1;
				break;
			case 16: Rho = temp_numcast; break;
			case 17: Rho1 = temp_numcast; break;
			case 18: Rho2 = temp_numcast; break;
			}
			break;
		default:
			throw std::runtime_error("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nUnknow input option, please use:\n --m1s --m2s --m3s --m4s --q10 --q20 --q21 --q30 --q31 --q32 --rho --rho1 --rho2\nOR\n--m1s --m2s --m3s --m4s --p1s --p2s --p3s --p4s --p12s --p23s --rho --rho1 --rho2\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
		}
	}
	if (loopToolInput_Indicator == 1 && xloopsInput_Indicator == 1)
		throw std::runtime_error("Please use only one type of input (xloop (q) /looptool (p))");
	else if(loopToolInput_Indicator == 1){
		ex p12, p23,p13;
		p12 = (p12s - p1s - p2s)/2;
		p23 = (p23s - p2s - p3s)/2;
		p13 = (p2s + p4s - p12s - p23s)/2;
		
		q10 = sqrt(p1s);
		q20 = p12/sqrt(p1s) + sqrt(p1s);
		q21 = sqrt( pow(p12/sqrt(p1s) + sqrt(p1s),2) - p12s );
		q30 = p13/sqrt(p1s) + p12/sqrt(p1s) + sqrt(p1s);
		q31 = (q30*q20 - p12s - p13 - p23)
			/
			q21;
		q32 = sqrt( pow(q30,2) - pow(q31,2) - p4s );
		
		/*input m2s, m3s, m4s, m1s for xloop*/
		m = lst(m2, m3, m4, m1);
	}
	else{
		m = lst(m1, m2, m3, m4);
	}
	q = lst(q10, q20, q21, q30, q31, q32);
}
