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
**	Level 7 functions
**
********************************************************************************
**
** Historial Log:
**	Date		Version	Author		Description
**	_________	_______	_________	________________________________
**	20090214	1.0	Quan Phan	Create this file
*******************************************************************************/

#include "ginac/ginac.h"
#include "trm.h"
#include "D0.h"
#include "OneLoop4Pt.h"

using namespace GiNaC;
namespace xloops{
	ex fn_1Loop4Pt(){
		ex return_value;
		mat_msquare[0]=symbol("m1");mat_msquare[1]=symbol("m2");mat_msquare[2]=symbol("m3");mat_msquare[3]=symbol("m4");
		mat_q[0][0]=symbol("q10");mat_q[0][1]=0;mat_q[0][2]=0;mat_q[0][3]=0;
		mat_q[1][0]=symbol("q20");mat_q[1][1]=symbol("q21");mat_q[1][2]=0;mat_q[1][3]=0;
		mat_q[2][0]=symbol("q30");mat_q[2][1]=symbol("q31");mat_q[2][2]=symbol("q32");mat_q[2][3]=0;
		
		trm_init();
		return_value = fn_D0();
		return return_value;
	}

	ex fn_1Loop4Pt(const ex& q_, const ex& m_, const ex& Rho_, const ex& Rho1_, const ex& Rho2_){
		ex return_value;
		mat_msquare[0]=m_.op(0);	mat_msquare[1]=m_.op(1);	mat_msquare[2]=m_.op(2);	mat_msquare[3]=m_.op(3);
		mat_q[0][0]=q_.op(0);mat_q[0][1]=0;mat_q[0][2]=0;mat_q[0][3]=0;
		mat_q[1][0]=q_.op(1);mat_q[1][1]=q_.op(2);mat_q[1][2]=0;mat_q[1][3]=0;
		mat_q[2][0]=q_.op(3);mat_q[2][1]=q_.op(4);mat_q[2][2]=q_.op(5);mat_q[2][3]=0;
		input_Rho = Rho_; input_Rho1 = Rho1_; input_Rho2 = Rho2_;
			//Rho = Rho_; Rho1 = Rho1_; Rho2 = Rho2_;
		trm_init();
		return_value = fn_D0();
		return return_value;
	}
}// Namespace xloops
