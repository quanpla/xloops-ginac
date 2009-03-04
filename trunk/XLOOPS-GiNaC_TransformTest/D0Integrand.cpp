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
**	Calculate the Integrands
**
********************************************************************************
**
** Historial Log:
**	Date		Version	Author		Description
**	_________	_______	_________	________________________________
*******************************************************************************/
//#define _AGK_debugmode

#include "ginac/ginac.h"

#include "trm.h"
#include "my_fns.h"
#include "trm2F.h"

#include "lev1.h"
#include "lev2.h"
#include "lev3.h"
#include "D0Integrand.h"

using namespace GiNaC;

namespace xloops{
	ex D0_integrand1(){
		xloopsGiNaC_calc_lev1();
		// integrand1 = 1/(P1*P2*P3*P4)
		ex l0 = realsymbol("l0"), l1 = realsymbol("l1"), l2 = realsymbol("l2"), lortho = realsymbol("lortho");
	
		ex P1, P2, P3, P4;
		ex integrand1;
		
		P1 = pow(l0 + mat_q[0][0], 2) - pow(l1, 2) - pow(l2, 2) - pow(lortho, 2) - mat_msquare[0] + I*Rho;
		P2 = pow(l0 + mat_q[1][0], 2) - pow(l1 + mat_q[1][1], 2) - pow(l2, 2) - pow(lortho, 2) - mat_msquare[1] + I*Rho;
		P3 = pow(l0 + mat_q[2][0], 2) - pow(l1 + mat_q[2][1], 2) - pow(l2 + mat_q[2][2], 2) - pow(lortho, 2) - mat_msquare[2] + I*Rho;
		P4 = pow(l0, 2) - pow(l1, 2) - pow(l2, 2) - pow(lortho, 2) - mat_msquare[3] + I*Rho;
		
		integrand1 = 1.0 / (P1*P2*P3*P4);
		
		return integrand1;
	}
	
	ex D0_integrand9(){
		xloopsGiNaC_calc_lev1();
		
		// integrand9 = 1/(l0^2 - l1^2 -...) / PI(alk.l0 + blk.l1 + ...)
		ex integrand9 = 0;
		
		ex l0 = realsymbol("l0"), l1 = realsymbol("l1"), l2 = realsymbol("l2"), lortho = realsymbol("lortho");
		
		int l, k;
		
		for(k = 0; k<4; k++){
			ex term1 = (pow(l0, 2) - pow(l1, 2) - pow(l2, 2) - pow(lortho, 2) - mat_msquare[k] + I*Rho);
			ex term2 = 1;
			for(l = 0; l<4; l++){
				if(l!=k){
					term2 *=(mat_a[l][k]*l0 + mat_b[l][k]*l1 + mat_c[l][k]*l2 + mat_d[l][k]);
				}
			}
			integrand9 += (1.0/term1) * (1.0/term2);
		}
		
		return integrand9;
	}
	
	ex D0_integrand12(){
		xloopsGiNaC_calc_lev1();
		xloopsGiNaC_calc_lev2(); // AC is in level 2 terms
		
		ex integrand12 = 0;
		
		ex x = realsymbol("x"), y = realsymbol("y"), z = realsymbol("z"), t = realsymbol("t");
		
		int l, k;
		
		for(k = 0; k<4; k++){
			ex term1 = 2.0*x*z - pow(z, 2) - pow(y, 2) - pow(t, 2) - mat_msquare[k] + I*Rho;
			
			ex term2 = 1;
			for(l = 0; l<4; l++){
				if(l!=k){
					term2 *=(mat_a[l][k]*z + mat_b[l][k]*y + mat_AC[l][k]*x + mat_d[l][k]);
				}
			}
			
			integrand12 += (1.0/term1) * (1.0/term2);
		}
		
		return integrand12;
	}
	
	ex D0_integrand1801(){ // equation 18, D0+ (01)
		xloopsGiNaC_calc_lev1();
		xloopsGiNaC_calc_lev2(); // AC_lk in lev 2 calculation
		
		ex integrand1801 = 0;
		
		ex x = realsymbol("x"), y = realsymbol("y"), z = realsymbol("z"), t = realsymbol("t");
		
		int m, l, k;
		
		for(k = 0; k<4; k++) for(l = 0; l<4; l++){
			if(l!=k){
				ex term1 = 1.0/mat_AC[l][k];
				ex term2 = 1;
				for(m=0; m<4; m++){
					if(m!=l && m!=k){
						term2 *= 1.0/(mat_A[m][l][k]*z + mat_B[m][l][k]*y + mat_C[m][l][k]);
						term2 *= mat_f[l][k] * (1.0 - myfn_delta(mat_AC[l][k]));
						term2 *= 1.0/( (1.0-2.0*mat_a[l][k]/mat_AC[l][k])*pow(z,2)
							     - 2.0*mat_b[l][k]*y*z/mat_AC[l][k]
							     - 2.0*mat_d[l][k]*z/mat_AC[l][k]
							     -pow(y, 2)
							     -pow(t, 2)
							     -mat_msquare[k]
							     + I*Rho);
					}
				}
				integrand1801 += term1 * term2;
			}
		}
		
		// because we get the integral for t limit = -infty..infty
		//	we will drop the *2 factor
		integrand1801 *= Pi*I;
		return integrand1801;
	}
	
	ex D0_integrand1802(){ // equation 18, D0-
		xloopsGiNaC_calc_lev1();
		xloopsGiNaC_calc_lev2(); // AC_lk in lev 2 calculation
		
		ex integrand1802 = 0;
		
		ex x = realsymbol("x"), y = realsymbol("y"), z = realsymbol("z"), t = realsymbol("t");
		
		int m, l, k;
		
		for(k = 0; k<4; k++) for(l = 0; l<4; l++){
			if(l!=k){
				ex term1 = 1.0/mat_AC[l][k];
				ex term2 = 1;
				for(m=0; m<4; m++){
					if(m!=l && m!=k){
						term2 *= 1.0/(mat_A[m][l][k]*z + mat_B[m][l][k]*y + mat_C[m][l][k]);
						/*the only different between D0+ and - is mat_f and mat_fminus*/
						term2 *= mat_fminus[l][k] * (1.0 - myfn_delta(mat_AC[l][k]));
						term2 *= 1.0/( (1.0-2.0*mat_a[l][k]/mat_AC[l][k])*pow(z,2)
								- 2.0*mat_b[l][k]*y*z/mat_AC[l][k]
								- 2.0*mat_d[l][k]*z/mat_AC[l][k]
								-pow(y, 2)
								-pow(t, 2)
								-mat_msquare[k]
								+ I*Rho);
					}
				}
				integrand1802 += term1 * term2;
			}
		}
		
		// the integrate is for limit of t = -infty..infty => no *2 factor
		integrand1802 *= - Pi*I;
		return integrand1802;
	}
	
	
	ex D0_integrand(int eqnum, const ex& q_, const ex& m_, const ex& Rho_, const ex& Rho1_, const ex& Rho2_){
		ex return_value;
		
		// init variables
		mat_msquare[0]=m_.op(0); mat_msquare[1]=m_.op(1); mat_msquare[2]=m_.op(2); mat_msquare[3]=m_.op(3);
		mat_q[0][0]=q_.op(0);mat_q[0][1]=0;mat_q[0][2]=0;mat_q[0][3]=0;
		mat_q[1][0]=q_.op(1);mat_q[1][1]=q_.op(2);mat_q[1][2]=0;mat_q[1][3]=0;
		mat_q[2][0]=q_.op(3);mat_q[2][1]=q_.op(4);mat_q[2][2]=q_.op(5);mat_q[2][3]=0;
		input_Rho = Rho_; input_Rho1 = Rho1_; input_Rho2 = Rho2_;
		Rho = Rho_; Rho1 = Rho1_; Rho2 = Rho2_;
		
		// return the desire integrand;
		switch(eqnum){
			case 1: return_value = D0_integrand1(); break;
			case 9: return_value = D0_integrand9(); break;
			case 12: return_value = D0_integrand12(); break;
			case 1801/*D0+*/: return_value = D0_integrand1801(); break;
			case 1802/*D0-*/: return_value = D0_integrand1802(); break;			
			default: return_value = 0;
		}
		
		return return_value;
	}
}// Namespace xloops
