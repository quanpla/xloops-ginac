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
#include "extrm.h"
#include "my_fns.h"
#include "trmchk.h"
#include "lev1.h"
#include "lev2.h"
#include "lev3.h"
#include "lev4.h"
#include "lev5.h"

using namespace GiNaC;
namespace xloops{

	ex fn_OPlus (int n, int m, int l, int k){
		ex OPlus_nmlk;
		//init
		ex f = fn_f(l, k), fminus = fn_fminus(l, k), g = fn_g(m, l, k), gminus = fn_gminus(m, l, k),
			F = fn_F(n, m, l, k), beta = fn_beta(m, l, k), P = fn_P(m, l, k), z1beta = fn_z1beta(m, l, k), z2beta = fn_z2beta(m, l, k), z1phi = fn_z1phi(m, l, k), z2phi = fn_z2phi(m, l, k),
		phi = fn_phi(m, l, k);
		check0denom(beta, "OPlus", n, m,l ,k);

	 // calculate
		OPlus_nmlk = (f*g + fminus*g)*log(F/beta);
		OPlus_nmlk += -2.0*Pi*I*( f*g + fminus*g )*my_step(imag_part(-P*z1beta/beta))*my_step(imag_part(z2beta));
		OPlus_nmlk += 2.0*Pi*I*f*g*my_step(imag_part(-P*phi*z1phi))*my_step(imag_part(z2phi));
		OPlus_nmlk += -2.0*Pi*I*f*gminus*my_step(imag_part(-P*phi*z1phi))*my_step(-imag_part(z2phi));
	
		return OPlus_nmlk;
	}

	ex fn_OMinus (int n, int m, int l, int k){
		ex OMinus_nmlk;
		//init
		ex f = fn_f(l, k), fminus = fn_fminus(l, k), g = fn_g(m, l, k), gminus = fn_gminus(m, l, k),
			F = fn_F(n, m, l, k), beta = fn_beta(m, l, k), P = fn_P(m, l, k), z1beta = fn_z1beta(m, l, k), z2beta = fn_z2beta(m, l, k), z1phi = fn_z1phi(m, l, k), z2phi = fn_z2phi(m, l, k),
		phi = fn_phi(m, l, k);
		check0denom(beta, "OPlus", n, m,l ,k);

	 // calculate
		OMinus_nmlk = -fminus*gminus*log(F/beta) - f*gminus*log(-F/beta);
		OMinus_nmlk += 2.0*Pi*I*fminus*gminus*my_step(imag_part(-P*z1beta/beta))*my_step(imag_part(z2beta));
		OMinus_nmlk += -2.0*Pi*I*f*gminus*my_step(imag_part(-P*z1beta/beta))*my_step(-imag_part(z2beta));
		OMinus_nmlk += -2.0*Pi*I*(fminus*gminus + fminus*g)*my_step(imag_part(-P*phi*z1phi))*my_step(imag_part(z2phi));
	
		return OMinus_nmlk;
	}

}// Namespace xloops
