/*******************************************************************************
**
**      xloop-ginacs Project
**      decompose log
**
**
**      HCMUNS,
**      Mar 200
**
**
**      Author(s):      Son, D. H.      (sondo01@gmail.com)
**                      Khiem, Phan     (phanhongkhiem@gmail.com)
**                      Quan, Phan      (anhquan.phanle@gmail.com)
********************************************************************************
**
** Historial Log:
**      Date            Version Author          Description
**      _________       _______ _________       ________________________________
*******************************************************************************/
#ifndef __XLOOPS_ONELOOP_4PT_LOGDECOMPOSE_H__
#define __XLOOPS_ONELOOP_4PT_LOGDECOMPOSE_H__


#include "ginac/ginac.h"

using namespace GiNaC;

namespace xloops{
ex log_NPoint_82_plus(int k, int l, int m, int n, const ex &sigma, const ex &z);
ex log_NPoint_82_minus(int k, int l, int m, int n, const ex &sigma, const ex &z);
ex log_Khiem_46_plus(int k, int l, int m, int n, const ex &sigma, const ex &z);
ex log_Khiem_46_minus(int k, int l, int m, int n, const ex &sigma, const ex &z);
ex log_original_plus(int k, int l, int m, int n, const ex &sigma, const ex &z);
ex log_original_minus(int k, int l, int m, int n, const ex &sigma, const ex &z);
ex fn_S(int k, int l, int m, int n, const ex &sigma, const ex &z);
ex fn_z1sigma(int k, int l, int m, int n, const ex &sigma);
ex fn_z2sigma(int k, int l, int m, int n, const ex &sigma);
} // namespace xloops

#endif  // __XLOOPS_ONELOOP_4PT_LOGDECOMPOSE_H__
