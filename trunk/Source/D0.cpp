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
//#define _AGK_debugmode

#include "ginac/ginac.h"

#include "extrm.h"
#include "my_fns.h"
#include "trmchk.h"

#include "lev1.h"
#include "lev2.h"
#include "lev3.h"
#include "lev4.h"
#include "lev5.h"
#include "lev6.h"
#include "lev7.h"

#include "D0.h"

using namespace GiNaC;

namespace xloops{
	void trm_init(){
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
					mat_z1phi[k][l][m] = 0;
					mat_z2phi[k][l][m] = 0;
					mat_z3phi[k][l][m] = 0;
					mat_z4phi[k][l][m] = 0;
					mat_z1beta[k][l][m] = 0;
					mat_z2beta[k][l][m] = 0;
					mat_z3beta[k][l][m] = 0;
					mat_z4beta[k][l][m] = 0;
					for(n=0;n<4;n++){
						mat_F[k][l][m][n] = 0; mat_F_im[k][l][m][n] = 0;
						mat_T1[k][l][m][n] = 0;
						mat_T2[k][l][m][n] = 0;
						mat_T3[k][l][m][n] = 0;
						mat_T4[k][l][m][n] = 0;
						mat_OPlus[k][l][m][n] = 0;
						mat_OMinus[k][l][m][n] = 0;
					}
				}
			}
		}
		///*Calculate Level 1 Terms*/
#ifdef _AGK_debugmode
		printf("Level 1 Calculating. Please Wait a Long Moment ...\n");
#endif
		for(k = 0; k<4; k++) for(l = 0; l<4; l++) if (l!=k){
#ifdef _AGK_debugmode
			printf("%d-%d\n", k+1, l+1);
#endif
			mat_a[l][k] = fn_a(l, k);
			mat_b[l][k] = fn_b(l, k);
			mat_c[l][k] = fn_c(l, k);
			mat_d[l][k] = fn_d(l, k); mat_d_im[l][k] = fn_d_im(l, k); mat_d_re[l][k] = fn_d_re(l, k); mat_d_conj[l][k] = fn_d_conj(l, k);
		}
		///*Calculate Level 2 Terms*/
	
#ifdef _AGK_debugmode
		printf("Level 2 Calculating. Please Wait a Long Moment ...\n");
#endif
	
		for(k = 0; k<4; k++) for(l = 0; l<4; l++) if (l!=k){
#ifdef _AGK_debugmode
			printf("%d-%d\n", k+1, l+1);
#endif
			mat_AC[l][k] = fn_AC(l, k);
			mat_alpha[l][k] = fn_alpha(l, k);
		}
		///*Calculate Level 3 Terms*/
	
#ifdef _AGK_debugmode
		printf("Level 3 Calculating. Please Wait a Long Moment ...\n");
#endif
	
		for(k = 0; k<4; k++) for(l = 0; l<4; l++) if (l!=k){
#ifdef _AGK_debugmode
			printf("%d-%d\n", k+1, l+1);
#endif
			mat_f[l][k] = fn_f(l, k);
			mat_fminus[l][k] = fn_fminus(l,k);
		}
		for(k = 0; k<4; k++) for(l = 0; l<4; l++) if (l!=k) for(m = 0; m<4; m++) if(m!=l && m!=k){
#ifdef _AGK_debugmode
			printf("%d-%d-%d\n", k+1, l+1, m+1);
#endif
			mat_A[m][l][k] = fn_A(m, l, k);
			mat_B[m][l][k] = fn_B(m, l, k);
			mat_C[m][l][k] = fn_C(m, l, k); mat_C_im[m][l][k] = fn_C_im(m, l, k); mat_C_re[m][l][k] = fn_C_re(m, l, k); mat_C_conj[m][l][k] = fn_C_conj(m, l, k);
			mat_D[m][l][k] = fn_D(m, l, k);
		}
		///*Calculate Level 4 Terms*/
	
#ifdef _AGK_debugmode
		printf("Level 4 Calculating. Please Wait a Long Moment ...\n");
#endif
	
		for(k = 0; k<4; k++) for(l = 0; l<4; l++) if (l!=k) for(m = 0; m<4; m++) if(m!=l && m!=k){
#ifdef _AGK_debugmode
			printf("%d-%d-%d\n", k+1, l+1, m+1);
#endif
			mat_beta[m][l][k]	=fn_beta(m, l, k);
			mat_phi[m][l][k]	=fn_phi(m, l, k);
			mat_g[m][l][k]		=fn_g(m, l, k);
			mat_gminus[m][l][k]	=fn_gminus(m, l, k);
		}
		///*Calculate Level 5 Terms*/
	
#ifdef _AGK_debugmode
		printf("Level 5 Calculating. Please Wait a Long Moment ...\n");
#endif
		for(k=0; k<4; k++) for(l=0; l<4; l++) if (l!=k) for(m = 0; m<4; m++) if(m!=l && m!=k){
#ifdef _AGK_debugmode
			printf("%d-%d-%d\n", k+1, l+1, m+1);
#endif
			mat_Q[m][l][k]	=fn_Q(m, l, k); mat_Q_im[m][l][k]=fn_Q_im(m, l, k); mat_Q_re[m][l][k]=fn_Q_re(m, l, k); mat_Q_conj[m][l][k]=fn_Q_conj(m, l, k);
			mat_P[m][l][k]	=fn_P(m, l, k);
			mat_E[m][l][k]	=fn_E(m, l, k); mat_E_im[m][l][k]=fn_E_im(m, l, k); mat_E_re[m][l][k]=fn_E_re(m, l, k);
		}
		for(k=0; k<4; k++) for(l=0; l<4; l++) if (l!=k) for(m=0; m<4; m++) if(m!=l && m!=k) for(n=0; n<4; n++) if(n!=m && n!=l && n!=k){
#ifdef _AGK_debugmode
			printf("%d-%d-%d-%d\n", k+1, l+1, m+1, n+1);
#endif
			mat_F[n][m][l][k] =fn_F(n, m, l, k); mat_F_im[n][m][l][k] =fn_F_im(n, m, l, k);
		}
	
		///*Calculate Level 6 Terms*/
#ifdef _AGK_debugmode
		printf("Level 6 Calculating. Please Wait a Long Moment ...\n");
#endif
		for(k = 0; k<4; k++) for(l = 0; l<4; l++) if (l!=k) for(m = 0; m<4; m++) if(m!=l && m!=k){
#ifdef _AGK_debugmode
			printf("%d-%d-%d\n", k+1, l+1, m+1);
#endif
			mat_A0[m][l][k]	=fn_A0(m, l, k);
			mat_B0[m][l][k]	=fn_B0(m, l, k);
			mat_C0[m][l][k]	=fn_C0(m, l, k);
			mat_z1phi[m][l][k]	=fn_z1phi(m, l, k);
			mat_z2phi[m][l][k]	=fn_z2phi(m, l, k);
			mat_z3phi[m][l][k]	=fn_z3phi(m, l, k);
			mat_z4phi[m][l][k]	=fn_z4phi(m, l, k);
			mat_z1beta[m][l][k]	=fn_z1beta(m, l, k);
			mat_z2beta[m][l][k]	=fn_z2beta(m, l, k);
			mat_z3beta[m][l][k]	=fn_z3beta(m, l, k);
			mat_z4beta[m][l][k]	=fn_z4beta(m, l, k);
		
		}
		for(k = 0; k<4; k++) for(l = 0; l<4; l++) if (l!=k) for(m = 0; m<4; m++) if(m!=l && m!=k) for(n = 0; n<4; n++) if(n!=m && n!=l && n!=k){
#ifdef _AGK_debugmode
			printf("%d-%d-%d-%d\n", k+1, l+1, m+1, n+1);
#endif
			mat_T1[n][m][l][k] = fn_T1(n, m, l, k);
			mat_T2[n][m][l][k] = fn_T2(n, m, l, k);
			mat_T3[n][m][l][k] = fn_T3(n, m, l, k);
			mat_T4[n][m][l][k] = fn_T4(n, m, l, k);
		}
	
		///*Calculate Level 7 Terms*/
#ifdef _AGK_debugmode
		printf("Level 7 Calculating. Please Wait a Long Moment ...\n");
#endif
		for(k = 0; k<4; k++) for(l = 0; l<4; l++) if (l!=k) for(m = 0; m<4; m++) if(m!=l && m!=k) for(n = 0; n<4; n++) if(n!=m && n!=l && n!=k){
#ifdef _AGK_debugmode
			printf("%d-%d-%d-%d\n", k+1, l+1, m+1, n+1);
#endif
			mat_OPlus[n][m][l][k] = fn_OPlus(n, m, l, k);
			mat_OMinus[n][m][l][k] = fn_OMinus(n, m, l, k);
		}
	}
	
	
	ex fn_D0(){
		/** **********************************************************************
		 ** This is the main integration, return the D0 of calculation.
		 *************************************************************************/
/// Calculate Level 8 terms
		int k, l, m, n;
#ifdef _AGK_debugmode
		printf("Level 8 Calculating. Please Wait a Long Moment ...\n");
#endif
		D0 = 0;
		
		for(k = 0; k<4; k++) for(l = 0; l<4; l++) if (l!=k) for(m = 0; m<4; m++) if(m!=l && m!=k) for(n = 0; n<4; n++) if(n!=m && n!=l && n!=k){
#ifdef _AGK_debugmode
			printf("%d-%d-%d-%d\n", k+1, l+1, m+1, n+1);
#endif
			ex term1 = 0, term2 = 0;
			term1 = mat_B[m][l][k]/(mat_B[m][l][k]*mat_A[n][l][k] - mat_B[n][l][k]*mat_A[m][l][k]);
			term1 *=(1.0-myfn_delta(mat_AC[l][k]))*(1.0-myfn_delta(mat_B[m][l][k]))*abs(1.0-mat_beta[m][l][k]*mat_phi[m][l][k])*(-1.0/mat_P[m][l][k]);
		
			term2 = /*1*/mat_OPlus[n][m][l][k]*fn_R(-mat_T1[n][m][l][k],-mat_T2[n][m][l][k]);
			term2 += /*2*/-mat_f[l][k]*mat_g[m][l][k]*fn_LogAG((1.0-mat_beta[m][l][k]*mat_phi[m][l][k])/mat_beta[m][l][k], (mat_F[n][m][l][k]+I*Rho*mat_beta[m][l][k])/mat_beta[m][l][k], -mat_T1[n][m][l][k], -mat_T2[n][m][l][k]);
			term2 += /*3*/-mat_f[l][k]*mat_gminus[m][l][k]*fn_LogAG((1.0-mat_beta[m][l][k]*mat_phi[m][l][k])/(-mat_beta[m][l][k]), (mat_F[n][m][l][k]+I*Rho*mat_beta[m][l][k])/(-mat_beta[m][l][k]), -mat_T1[n][m][l][k], -mat_T2[n][m][l][k]);
			term2 += /*4*/-(mat_f[l][k]*mat_gminus[m][l][k] + mat_f[l][k]*mat_g[m][l][k])*fn_LogAG(-1.0/mat_beta[m][l][k], mat_z1beta[m][l][k]/mat_beta[m][l][k], -mat_T1[n][m][l][k], -mat_T2[n][m][l][k]);
			term2 += /*5*/-(mat_f[l][k]*mat_gminus[m][l][k] + mat_f[l][k]*mat_g[m][l][k])*fn_LogAG(1.0, -mat_z2beta[m][l][k], -mat_T1[n][m][l][k], -mat_T2[n][m][l][k]);
			term2 += /*6*/mat_f[l][k]*mat_g[m][l][k]*fn_LogAG(-mat_phi[m][l][k], mat_phi[m][l][k]*mat_z1phi[m][l][k], -mat_T1[n][m][l][k], -mat_T2[n][m][l][k]);
			term2 += /*7*/mat_f[l][k]*mat_g[m][l][k]*fn_LogAG(1.0, -mat_z2phi[m][l][k], -mat_T1[n][m][l][k], -mat_T2[n][m][l][k]);
			term2 += /*8*/mat_f[l][k]*mat_gminus[m][l][k]*fn_LogAG(mat_phi[m][l][k], -mat_phi[m][l][k]*mat_z1phi[m][l][k], -mat_T1[n][m][l][k], -mat_T2[n][m][l][k]);
			term2 += /*9*/mat_f[l][k]*mat_gminus[m][l][k]*fn_LogAG(1.0, -mat_z2phi[m][l][k], -mat_T1[n][m][l][k], -mat_T2[n][m][l][k]);
			term2 += /*10*/mat_fminus[l][k]*mat_g[m][l][k]*fn_LogAG(1.0, (mat_Q[m][l][k] + I*Rho)/mat_P[m][l][k], -mat_T1[n][m][l][k], -mat_T2[n][m][l][k]);
			term2 += /*11*/-mat_f[l][k]*mat_gminus[m][l][k]*fn_LogAG(1.0, (mat_Q[m][l][k] - I*Rho)/mat_P[m][l][k], -mat_T1[n][m][l][k], -mat_T2[n][m][l][k]);
		
			term2 = /*12*/mat_OMinus[n][m][l][k]*fn_R(-mat_T3[n][m][l][k],-mat_T4[n][m][l][k]);
			term2 += /*13*/(mat_fminus[l][k]*mat_gminus[m][l][k] + mat_fminus[l][k]*mat_g[m][l][k])*fn_LogAG(-(1.0-mat_beta[m][l][k]*mat_phi[m][l][k])/mat_beta[m][l][k], mat_F[n][m][l][k]/mat_beta[m][l][k], -mat_T3[n][m][l][k], -mat_T4[n][m][l][k]);
			term2 += /*14*/mat_fminus[l][k]*mat_gminus[m][l][k]*fn_LogAG(-1.0/mat_beta[m][l][k], mat_z3beta[m][l][k]/mat_beta[m][l][k], -mat_T3[n][m][l][k], -mat_T4[n][m][l][k]);
			term2 += /*15*/mat_fminus[l][k]*mat_gminus[m][l][k]*fn_LogAG(1.0, mat_z4beta[m][l][k], -mat_T3[n][m][l][k], -mat_T4[n][m][l][k]);
			term2 += /*16*/mat_f[l][k]*mat_gminus[m][l][k]*fn_LogAG(1.0/mat_beta[m][l][k],-mat_z3beta[m][l][k]/mat_beta[m][l][k], -mat_T3[n][m][l][k], -mat_T4[n][m][l][k]);
			term2 += /*17*/mat_f[l][k]*mat_gminus[m][l][k]*fn_LogAG(1.0, mat_z4beta[m][l][k], -mat_T3[n][m][l][k], -mat_T4[n][m][l][k]);
			term2 += /*18*/-(mat_fminus[l][k]*mat_gminus[m][l][k] + mat_f[l][k]*mat_gminus[m][l][k])*fn_LogAG(-mat_phi[m][l][k], mat_phi[m][l][k]*mat_z3phi[m][l][k], -mat_T3[n][m][l][k], -mat_T4[n][m][l][k]);
			term2 += /*19*/-(mat_fminus[l][k]*mat_gminus[m][l][k] + mat_f[l][k]*mat_gminus[m][l][k])*fn_LogAG(1.0, -mat_z4phi[m][l][k], -mat_T3[n][m][l][k], -mat_T4[n][m][l][k]);
			term2 += /*20*/mat_fminus[l][k]*mat_g[m][l][k]*fn_LogAG(-1.0, (mat_Q[m][l][k] + I*Rho)/mat_P[m][l][k], -mat_T3[n][m][l][k], -mat_T4[n][m][l][k]);
			term2 += /*21*/-mat_f[l][k]*mat_gminus[m][l][k]*fn_LogAG(-1.0, (mat_Q[m][l][k] - I*Rho)/mat_P[m][l][k], -mat_T3[n][m][l][k], -mat_T4[n][m][l][k]);
			
			D0 += term1 * term2;
		}
		printf("Finish calculating. Multiply with I.Pi^2.\n");
		D0 *= I*pow(Pi,2);
		//printf("Evaluating.\n");
		//D0 = my_evalf(D0);
#ifdef _AGK_debugTestCase1
		exportTerms2File("allLev.txt");
#endif
		return D0;
	}

}// Namespace xloops
