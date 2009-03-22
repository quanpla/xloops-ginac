/*******************************************************************************
**
** xloop-ginacs Project
** 1Loop4Pt implementation
**
**
** HCMUNS,
** June 2008
**
**
** Author(s): Son, D. H. (sondo01@gmail.com)
** Khiem, Phan (phanhongkhiem@gmail.com)
** Quan, Phan (anhquan.phanle@gmail.com)
********************************************************************************
**
** This file contains all the terms necessary to calculate D0 in Npoint.
** Of course it must
**
********************************************************************************
**
** Historial Log:
** Date Version Author Description
** _________ _______ _________ ________________________________
**
*******************************************************************************/


#ifndef __XLOOPS_ONELOOP_4PT_LEV4_H__
#define __XLOOPS_ONELOOP_4PT_LEV4_H__

#include <ginac/ginac.h>

using namespace std;

namespace xloops{
	
	ex fn_z1phi (int m, int l, int k);
	ex fn_z2phi (int m, int l, int k);
	ex fn_z3phi (int m, int l, int k);
	ex fn_z4phi (int m, int l, int k);
	ex fn_z1beta (int m, int l, int k);
	ex fn_z2beta (int m, int l, int k);
	ex fn_z3beta (int m, int l, int k);
	ex fn_z4beta (int m, int l, int k);

	ex fn_T1 (int n, int m, int l, int k);
	ex fn_T2 (int n, int m, int l, int k);
	ex fn_T3 (int n, int m, int l, int k);
	ex fn_T4 (int n, int m, int l, int k);
	
	void xloopsGiNaC_calc_lev4(); // calculate all level 4 terms
} // Namespace xloops

#endif // __XLOOPS_ONELOOP_4PT_LEV4_H__
