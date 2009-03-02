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

#include "lev1.h"

#include "D0Integrand.h"

using namespace GiNaC;

namespace xloops{
	ex l0 = realsymbol("l0"), l1 = realsymbol("l1"), l2 = realsymbol("l2"), lortho = realsymbol("lortho");
	
	ex D0_integrand1(){
		init_lev1();
		// integrand1 = 1/(P1*P2*P3*P4)
		
		ex P1, P2, P3, P4;
		ex integrand1;
		
		P1 = pow(l0 + mat_q[0][0], 2) - pow(l1, 2) - pow(l2, 2) - pow(lortho, 2) - mat_msquare[0] + I*Rho;
		P2 = pow(l0 + mat_q[1][0], 2) - pow(l1 + mat_q[1][0], 2) - pow(l2, 2) - pow(lortho, 2) - mat_msquare[1] + I*Rho;
		P3 = pow(l0 + mat_q[2][0], 2) - pow(l1 + mat_q[2][1], 2) - pow(l2 + mat_q[2][2], 2) - pow(lortho, 2) - mat_msquare[2] + I*Rho;
		P4 = pow(l0, 2) - pow(l1, 2) - pow(l2, 2) - pow(lortho, 2) - mat_msquare[3] + I*Rho;
		
		integrand1 = 1.0 / (P1*P2*P3*P4);
		
		return integrand1;
	}
	
	ex D0_integrand9(){
		init_lev1();
	}
	
	ex D0_integrand12(){
		init_lev1();
	}
	
	ex D0_integrand18(){
		init_lev1();
	}
	
	
	ex D0_integrand(int eqnum, const ex& q_, const ex& m_, const ex& Rho_, const ex& Rho1_, const ex& Rho2_){
		ex return_value;
		
		// init variables
		mat_msquare[0]=m_.op(0);	mat_msquare[1]=m_.op(1);	mat_msquare[2]=m_.op(2);	mat_msquare[3]=m_.op(3);
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
			case 18: return_value = D0_integrand18(); break;
			default: return_value = 0;
		}
		
		return return_value;
	}
}// Namespace xloops
