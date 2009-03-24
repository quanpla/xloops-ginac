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

ex fn_delta89 (int m, int l, int k);
ex fn_X189 (int m, int l, int k);
ex fn_X289 (int m, int l, int k);
ex fn_etaplus90 (int m, int l, int k, const ex &z);
ex fn_etaminus90 (int m, int l, int k, const ex &z);

} // Namespace xloops

#endif // __XLOOPS_ONELOOP_4PT_LEV4_H__
