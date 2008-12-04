
/* This is program testxloop.cpp, to compare to onloop4pt.cpp, we write this program */
/*       author: Phan Hong Khiem  
/*		 Phan Le Anh Quan
/*		 Do Hoang Son.                                                       */
/*       email: phanhongkhiem@gmail.com                                              */
/** *****************************************************************************
**
**	Execution Plan
**			+-------------------------------+
**			|				|
**			|            Input P ,m
**			|				|
**			+-------------------------------+
**					|
**					V
**			+-------------------------------+
**			|				|
**			|       calculate q{ij}	        |
**			|                       	|
**			|				|
**			+-------------------------------+
**                                      |
**                                      V
**			+-------------------------------+
**			|	    step1		|
**			| 				|
**			|      a_{lk},b_{l}{k},c_{lk}   |
**			|				|
**			+-------------------------------+
**                                      |
**                                      V
**			+-------------------------------+
**			|	   step2		|
**			|     A_{mlk}, B_{mlk},C_{mlk}  |
			|     ........F_{nmlk}......    |
**			|				|
**			+-------------------------------+
**                                      |
**                                      V
**			+-------------------------------+
**			|	   step 3		|
**			| 	call vegas calculate	|
**			|           D_{0}		|
**			|				|
**			+-------------------------------+
**
*******************************************************************************/
#include <iostream>
#include <getopt.h>
#include <sstream>
#include "ginac/ginac.h"
#include "ginac/flags.h"

using namespace std;
using namespace GiNaC;

// declare some function.
ex theta(ex );
ex delta(ex);
ex Etafunction (ex, ex);
ex R(ex, ex);
ex Gpositive (ex,ex);
ex Gnegative (ex,ex);
ex LnGpositive(ex, ex, ex);
ex LnGnegative(ex, ex, ex);

// input function
void input_option(int argc, char *argv[], lst &p, lst &m);

// main start point
int main(int argc, char *argv[])
{
	ex eta=1e-25;
	ex rho=2e-25;
// declare variable.
	int l,k,m,n;
	ex q10;
	ex q20,q21;
	ex q30,q31,q32;
	ex m1s,m2s,m3s,m4s;
	ex ms[5];
// declare a,b,c, d,AC,alpha  matrix.
	ex a[5][5];
	ex b[5][5];
	ex c[5][5];
	ex d[5][5];
	ex AC[5][5];
	ex alpha[5][5];
	ex f[5][5];
	ex fminus[5][5];
// declare A,B,C,D,beta,phi, Q,....f,g  matrix.
	ex A[5][5][5];
	ex B[5][5][5];
	ex C[5][5][5];
	ex D[5][5][5];
	ex beta[5][5][5];
	ex phi[5][5][5];
	ex Q[5][5][5];
	ex E[5][5][5];
	ex P[5][5][5];
	ex g[5][5][5];
	ex gminus[5][5][5];
	ex z1phi[5][5][5];
	ex z2phi[5][5][5];
	ex z1beta[5][5][5];
	ex z2beta[5][5][5];
	
// declare F matrix.
	ex F[5][5][5][5];
	ex Jacobian[5][5][5][5];
	ex T1[5][5][5][5];
	ex T2[5][5][5][5];
	
	
//cout<<"=================================== input p square and m square =============================================" << endl;
	// get the input
	lst m_, p_;
	input_option(argc, argv, p_, m_);
	
	ex p1s = p_.op(0);
	ex p2s = p_.op(1);
	ex p3s = p_.op(2);
	ex p4s = p_.op(3);
	ex p12s = p_.op(4);
	ex p23s = p_.op(5);
	m1s = m_.op(0);
	m2s = m_.op(1);
	m3s = m_.op(2);
	m4s = m_.op(3);
//cout<<"==================================== calculate q_{ij}---Lorent tranformation ================"<<endl;
	ex  p12 = (p12s - p1s - p2s)/2;
	ex  p23 = (p23s - p2s - p3s)/2;
	ex  p13 = (p2s + p4s - p12s - p23s)/2;
	
	q10 = sqrt(p1s);
	q20 = p12/sqrt(p1s) + sqrt(p1s);
	q21 = sqrt( pow(p12/sqrt(p1s) + sqrt(p1s),2) - p12s );
	q30 = p13/sqrt(p1s) + p12/sqrt(p1s) + sqrt(p1s);
	q31 = (q30*q20 - p12s - p13 - p23)/q21;	
	q32 = sqrt( pow(q30,2) - pow(q31,2) - p4s );
	
/* and   input q_{ij}*/
	
	ex q11=0, q12=0, q13=0, q22=0, q23=0, q33=0,q40=0, q41=0, q42=0, q43=0;
//cout<<"============================================================================================"<<endl;
//cout<<"=============================declare matrix for a{i}{j}====================================="<<endl;
	/* declare matrix for a{i}{j}*/
	a[1][1]=0.0;	
	a[1][2]=(q10-q20)*2.0;
	a[1][3]=(q10-q30)*2.0;
	a[1][4]=(q10-q40)*2.0;
	
	/*for l=2, k change*/
	
	a[2][1]=(q20-q10)*2.0;
	a[2][3]=(q20-q30)*2.0;
	a[2][4]=(q20-q40)*2.0;
	a[2][2]=0;
	/*for l=3, k change*/
	a[3][1]=(q30-q10)*2.0;
	a[3][2]=(q30-q20)*2.0;
	a[3][4]=(q30-q40)*2.0;
	a[3][3]=0;
	/*for l=4, k change*/
	a[4][1]=(q40-q10)*2.0;
	a[4][2]=(q40-q20)*2.0;
	a[4][3]=(q40-q30)*2.0;
	a[4][4]=0;
//cout<<"============================================================================================"<<endl;
//cout<<"=============================declare matrix for b{i}{j}====================================="<<endl;
	/*check b_{lk}=-2(q_{l1}-q_{k1})*/
	/*check for l=1,k change */
	b[1][1]=0;
	b[1][2]=-(q11-q21)*2.0;
	b[1][3]=-(q11-q31)*2.0;
	b[1][4]=-(q11-q41)*2.0;
	
	/*for l=2, k change*/
	
	b[2][1]=-(q21-q11)*2.0;
	b[2][2]=0.0;
	b[2][3]=-(q21-q31)*2.0;
	b[2][4]=-(q21-q41)*2.0;
	
	/*for l=3, k change */
	b[3][1]=-(q31-q11)*2.0;
	b[3][2]=-(q31-q21)*2.0;
	b[3][3]=0.0;
	b[3][4]=-(q31-q41)*2.0;
	/*for l=4,k change*/
	b[4][1]=-(q41-q11)*2.0;
	b[4][2]=-(q41-q21)*2.0;
	b[4][4]=0.0;
	b[4][3]=-(q41-q31)*2.0;
        /* cout b_{lk}*/
//cout<<"============================================================================================"<<endl;
//cout<<"=============================declare matrix for c{i}{j}====================================="<<endl;
	/*check c_{lk}=-2(ql2-qk2)*/
	/*for l=1 , k c[hange */
	c[1][1]=0;
	c[1][2]=-(q12-q22)*2.0;
	c[1][3]=-(q12-q32)*2.0;
	c[1][4]=-(q12-q42)*2.0;
	
	/*for l=2, k c[hange*/
	
	c[2][1]=-(q22-q12)*2.0;
	c[2][2]=0.0;
	c[2][3]=-(q22-q32)*2.0;
	c[2][4]=-(q22-q42)*2.0;
	
	/*for l=3, k c[hange */
	c[3][1]=-(q32-q12)*2.0;
	c[3][2]=-(q32-q22)*2.0;
	c[3][3]=0.0;
	c[3][4]=-(q32-q42)*2.0;
	/*for l=4,k c[hange*/
	c[4][1]=-(q42-q12)*2.0;
	c[4][2]=-(q42-q22)*2.0;
	c[4][4]=0.0;
	c[4][3]=-(q42-q32)*2.0;
//cout<<"==========================================================================================="<<endl;
//cout<<"==============================input d{lk}=================================================="<<endl;
	/*check d_{lk}=-2(ql2-qk2)*/
	/*for l=1 , k d[hange */
	d[1][1]=0;
	d[1][2]=pow(q10-q20,2)-pow(q21,2)-(m1s-m2s);
	d[1][3]=pow(q10-q30,2)-pow(q31,2)-pow(q32,2)-(m1s-m3s);
	d[1][4]=pow(q10,2)-(m1s-m4s);
	
	/*for l=2, k d[hange*/
	
	d[2][1]=pow(q20-q10,2)-pow(q21,2)-(m2s-m1s);
	d[2][2]=0.0;
	d[2][3]=pow(q20-q30,2) - pow(q21-q31,2) - pow(q32,2) - (m2s-m3s);
	d[2][4]=pow(q20,2)-pow(q21,2)-(m2s-m4s);
	
	/*for l=3, k d[hange */
	d[3][1]=pow(q30-q10,2)-pow(q31,2)-pow(q32,2)-(m3s-m1s);
	d[3][2]=pow(q30-q20,2)-pow(q31-q21,2)-pow(q32,2)-(m3s-m2s);
	d[3][3]=0.0;
	d[3][4]=pow(q30,2)-pow(q31,2)-pow(q32,2)-(m3s-m4s);
	/*for l=4,k d[hange*/
	d[4][1]=pow(q10,2)-(m4s-m1s);
	d[4][2]=pow(q20,2)-pow(q21,2)-(m4s-m2s);
	d[4][4]=0.0;
	d[4][3]=pow(q30,2)-pow(q31,2)-pow(q32,2)-(m4s-m3s);
	/* cout d_{lk}*/
//cout<<"============================================================================================="<<endl;
//cout<<"===================calculate AC_{lk}========================================================="<<endl;
	for(k = 1; k<=4; k++) for(l = 1; l<=4; l++) if(l!=k)
	{
		AC[k][l]=a[k][l]+c[k][l];
	} 
	
//cout<<"============================================================================================="<<endl;
//cout<<"========================calculate alpha_{lk}================================================="<<endl;
	for(k = 1; k<=4; k++) for(l = 1; l<=4; l++) if(l!=k)
	{
		if(AC[k][l]==0){
			cout<<"alpha_"<< l << k << "have denorminator is zero (:????"<<endl;
		}	
		else
		{		
			alpha[k][l]=b[k][l]/AC[k][l];
		}
	} 
	
//cout<<"==========================================================================================="<<endl;
//cout<<"========================calculate A_{mlk}=================================================="<<endl;
	
	for(k = 1; k<=4; k++) for(l = 1; l<=4; l++) if (l!=k) for(m = 1; m<=4; m++) if(m!=l && m!=k){
		if(AC[l][k]==0){
			cout<<"A_"<< m << l << k << "have denorminator is zero (:????"<<endl;
		}	
		else
		{		
			A[m][l][k]=a[m][k]-a[l][k]*AC[m][k]/AC[l][k];
		}
	}   
	
//cout<<"==========================================================================================="<<endl;
//cout<<"========================calculate B_{mlk}=================================================="<<endl;
	
	for(k = 1; k<=4; k++) for(l = 1; l<=4; l++) if (l!=k) for(m = 1; m<=4; m++) if(m!=l && m!=k){
		if(AC[l][k]==0){
			cout<<"B_"<< m << l << k << "have denorminator is zero (:????"<<endl;
		}	
		else
		{		
			B[m][l][k]=b[m][k]-b[l][k]*AC[m][k]/AC[l][k];
		}
	}   
	
	
	
//cout<<"==========================================================================================="<<endl;
//cout<<"========================calculate C_{mlk}=================================================="<<endl;
	
	for(k = 1; k<=4; k++) for(l = 1; l<=4; l++) if (l!=k) for(m = 1; m<=4; m++) if(m!=l && m!=k){
		if(AC[l][k]==0){
			cout<<"C_"<< m << l << k << "have denorminator is zero (:????"<<endl;
		}	
		else
		{		
			C[m][l][k]=d[m][k]-d[l][k]*AC[m][k]/AC[l][k];
		}
	}   
//cout<<"==========================================================================================="<<endl;
//cout<<"======================== calculate D_{mlk} ================================================"<<endl;
	for(k = 1; k<=4; k++) for(l = 1; l<=4; l++) if (l!=k) for(m = 1; m<=4; m++) if(m!=l && m!=k){
		if(AC[l][k]==0){
			cout<<"D_"<< m << l << k << "have denorminator is zero (:????"<<endl;
		}	
		else
		{
			D[m][l][k]=1.0-2.0*a[l][k]/AC[l][k]	+alpha[l][k]*alpha[l][k];
		}
	}
//cout<<"==========================================================================================="<<endl;
//cout<<"======================== calculate beta_{mlk} ================================================"<<endl;
	
	for(k = 1; k<=4; k++) for(l = 1; l<=4; l++) if (l!=k) for(m = 1; m<=4; m++) if(m!=l && m!=k){
		if(D[m][l][k]==0){
			cout<<"beta_"<< m << l << k << "have denorminator is zero (:????"<<endl;
		}	
		else
		{
			beta[m][l][k]=(
			                (A[m][l][k]/B[m][l][k]-alpha[l][k])
			                +sqrt((A[m][l][k]/B[m][l][k]-alpha[l][k])*(A[m][l][k]/B[m][l][k]-alpha[l][k])
			                      -D[m][l][k] +I*eta)
			              )/D[m][l][k];
		}
	}
//cout<<"==========================================================================================="<<endl;
//cout<<"======================== calculate phi_{mlk} ================================================"<<endl;
	
	for(k = 1; k<=4; k++) for(l = 1; l<=4; l++) if (l!=k) for(m = 1; m<=4; m++) if(m!=l && m!=k){
		
		phi[m][l][k]=
			(A[m][l][k]/B[m][l][k]-alpha[l][k])
			+sqrt((A[m][l][k]/B[m][l][k]-alpha[l][k])*(A[m][l][k]/B[m][l][k]-alpha[l][k])
			      -D[m][l][k] +I*eta);
		
		
	}
	
//cout<<"==========================================================================================="<<endl;
//cout<<"======================== calculate Q_{mlk} ================================================"<<endl;
	for(k = 1; k<=4; k++) for(l = 1; l<=4; l++) if (l!=k) for(m = 1; m<=4; m++) if(m!=l && m!=k){
		if(AC[l][k]==0 || B[m][l][k]==0){
			cout<<"Q_"<< m << l << k << "have denorminator is zero (:????"<<endl;
		}	
		else
		{
			Q[m][l][k]=-2.0*(C[m][l][k]/B[m][l][k]+d[l][k]*beta[m][l][k]/AC[l][k]);
		}
	}
//cout<<"==========================================================================================="<<endl;
//cout<<"======================== calculate E_{mlk} ================================================"<<endl;
	for(k = 1; k<=4; k++) for(l = 1; l<=4; l++) if (l!=k) for(m = 1; m<=4; m++) if(m!=l && m!=k){
		if(AC[l][k]==0 || B[m][l][k]==0){
			cout<<"E_"<< m << l << k << "have denorminator is zero (:????"<<endl;
		}	
		else
		{
			E[m][l][k]=-2.0*(C[m][l][k]*phi[m][l][k]/B[m][l][k]+d[l][k]/AC[l][k] );
		}
	}
//cout<<"==========================================================================================="<<endl;
//cout<<"======================== calculate P_{mlk} ================================================"<<endl;
	
	for(k = 1; k<=4; k++) for(l = 1; l<=4; l++) if (l!=k) for(m = 1; m<=4; m++) if(m!=l && m!=k){
		if( B[m][l][k]==0){
			cout<<"P_"<< m << l << k << "have denorminator is zero (:????"<<endl;
		}	
		else
		{
			P[m][l][k]=-2.0*( 
			                  (A[m][l][k]/B[m][l][k]-alpha[l][k])*(1+beta[m][l][k]*phi[m][l][k])
			                  -D[m][l][k]*beta[m][l][k]-phi[m][l][k]
			                );
		}
	}
	
//cout<<"==========================================================================================="<<endl;
//cout<<"======================== calculate f_{lk} ================================================"<<endl;
	for(k = 1; k<=4; k++) for(l = 1; l<=4; l++) if (l!=k){
		
		ex temp=imag_part(-d[l][k]/AC[l][k]);	
		
		if (temp>0) {
			f[l][k] = 2.0;
		}			
		if (temp==0) {
			f[l][k] =  1.0;
		}	
		if (temp <0 ){
			f[l][k] = 0.0;
		}
		
	}
//cout<<"==========================================================================================="<<endl;
//cout<<"======================== calculate fminus_{lk} ============================================"<<endl;
	for(k = 1; k<=4; k++) for(l = 1; l<=4; l++) if (l!=k)for(m = 1; m<=4; m++) if(m!=l && m!=k){
		
		ex temp=imag_part(-d[l][k]/AC[l][k]);	
		
		if (temp<0) {
			fminus[l][k] = 2.0;
		}			
		if (temp==0) {
			fminus[l][k] =  1.0;
		}	
		if (temp >0 ){
			fminus[l][k] = 0.0;
		}
		
	}
//cout<<"==========================================================================================="<<endl;
//cout<<"======================== calculate g_{lk} ================================================"<<endl;
	for(k = 1; k<=4; k++) for(l = 1; l<=4; l++) if (l!=k)for(m = 1; m<=4; m++) if(m!=l && m!=k){
		
		ex temp1=imag_part(-C[m][l][k]/B[m][l][k]);	
		
		if (temp1 > 0) {
			g[m][l][k] = 2.0;
		}			
		if (temp1==0) {
			g[m][l][k] =  1.0;
		}	
		if (temp1 <0 ){
			g[m][l][k] = 0.0;
		}
		
	}
//cout<<"==========================================================================================="<<endl;
//cout<<"======================== calculate gminus_{lk} ============================================"<<endl;
	
	for(k = 1; k<=4; k++) for(l = 1; l<=4; l++) if (l!=k)for(m = 1; m<=4; m++) if(m!=l && m!=k){
		
		ex temp1=imag_part(-C[m][l][k]/B[m][l][k]);	
		
		if (temp1<0) {
			gminus[m][l][k] = 2.0;
		}			
		if (temp1==0) {
			gminus[m][l][k] = 1.0;
		}	
		if (temp1 >0 ){
			gminus[m][l][k] = 0.0;
		}
		
	}
//cout<<"=========================calculate Z1phi and T2phi=================================================="<<endl
		for(k = 1; k<=4; k++) for(l = 1; l<=4; l++) if (l!=k) for(m = 1; m<=4; m++) if(m!=l && m!=k){
		   if(P[m][l][k]*phi[m][l][k]==0.0){
			cout<<"Z1phi_"<<m << l << k << "have denorminator is zero "<<endl;  }
		   else {
			z1phi[m][l][k]=((E[m][l][k]-Q[m][l][k]*phi[m][l][k])
			+sqrt( (E[m][l][k]-Q[m][l][k]*phi[m][l][k])*(E[m][l][k]-Q[m][l][k]*phi[m][l][k])-4.0*P[m][l][k]*phi[m][l][k]*(ms[k]-I*rho) ))/
				(2.0*P[m][l][k]*phi[m][l][k]);	  		  	 
		             }
			}
		 for(k = 1; k<=4; k++) for(l = 1; l<=4; l++) if (l!=k) for(m = 1; m<=4; m++) if(m!=l && m!=k)for(n= 1; n<=4; n++) if(n!=l && n!=k && n!=m){
    		  	 if(P[m][l][k]*phi[m][l][k]==0.0){
			cout<<"Z2phi_"<<m << l << k << "have denorminator is zero "<<endl;  }
		   else {
			z2phi[m][l][k]=((E[m][l][k]-Q[m][l][k]*phi[m][l][k])
			-sqrt( (E[m][l][k]-Q[m][l][k]*phi[m][l][k])*(E[m][l][k]-Q[m][l][k]*phi[m][l][k])-4.0*P[m][l][k]*phi[m][l][k]*(ms[k]-I*rho) ))/
				(2.0*P[m][l][k]*phi[m][l][k]);	  		  	 
		          }
			}	
//cout<<"=========================calculate Z1beta and T2beta=================================================="<<endl
		for(k = 1; k<=4; k++) for(l = 1; l<=4; l++) if (l!=k) for(m = 1; m<=4; m++) if(m!=l && m!=k){
		   if(P[m][l][k]*phi[m][l][k]==0.0){
			cout<<"Z1beta_"<<m << l << k << "have denorminator is zero "<<endl;  }
		   else {
			z1beta[m][l][k]=((E[m][l][k]-Q[m][l][k]/beta[m][l][k])
			+sqrt( (E[m][l][k]-Q[m][l][k]/beta[m][l][k])*(E[m][l][k]-Q[m][l][k]/beta[m][l][k])-4.0*(P[m][l][k]/beta[m][l][k])*(ms[k]-I*rho) ))/
				(2.0*P[m][l][k]/beta[m][l][k]);	  		  	 
		             }
			}
		 for(k = 1; k<=4; k++) for(l = 1; l<=4; l++) if (l!=k) for(m = 1; m<=4; m++) if(m!=l && m!=k)for(n= 1; n<=4; n++) if(n!=l && n!=k && n!=m){
    		  	 if(P[m][l][k]*phi[m][l][k]==0.0){
			cout<<"Z2beta_"<<m << l << k << "have denorminator is zero "<<endl;  }
		   else {
			z2beta[m][l][k]= ((E[m][l][k]-Q[m][l][k]/beta[m][l][k])
			-sqrt( (E[m][l][k]-Q[m][l][k]/beta[m][l][k])*(E[m][l][k]-Q[m][l][k]/beta[m][l][k])-4.0*(P[m][l][k]/beta[m][l][k])*(ms[k]-I*rho) ))/
				(2.0*P[m][l][k]/beta[m][l][k]);		  	 
		          }
			}	
	
	
//cout<<"==========================================================================================="<<endl;
//cout<<"======================== calculate E_{mlk} ================================================"<<endl;
	
	for(k = 1; k<=4; k++) for(l = 1; l<=4; l++) if (l!=k) for(m = 1; m<=4; m++) if(m!=l && m!=k)for(n= 1; n<=4; n++) if(n!=l && n!=k && n!=m){
		if(A[n][l][k]*B[m][l][k]-A[m][l][k]*B[n][l][k]==0){
			cout<<"F_"<<n<< m << l << k << "have denorminator is zero (:????"<<endl;
		}	
		else
		{
			F[n][m][l][k]=(C[n][l][k]*B[m][l][k]-C[m][l][k]*B[n][l][k])/(A[n][l][k]*B[m][l][k]-A[m][l][k]*B[n][l][k]);
		}
	}


//cout<<"=========================calculate T1 and T2=================================================="<<endl
		for(k = 1; k<=4; k++) for(l = 1; l<=4; l++) if (l!=k) for(m = 1; m<=4; m++) if(m!=l && m!=k)for(n= 1; n<=4; n++) if(n!=l && n!=k && n!=m){
    		  	 T1[n][m][l][k]=(Q[m][l][k]+P[m][l][k]*F[n][m][l][k]-beta[m][l][k]*E[m][l][k]
				+sqrt(pow(Q[m][l][k]+P[m][l][k]*F[n][m][l][k]-beta[m][l][k]*E[m][l][k],2.0)-4.0*P[m][l][k]*(Q[m][l][k]*F[n][m][l][k] 
						+beta[m][l][k]*ms[k]-I*beta[m][l][k]*rho )) )/(-P[m][l][k]);
			}

		for(k = 1; k<=4; k++) for(l = 1; l<=4; l++) if (l!=k) for(m = 1; m<=4; m++) if(m!=l && m!=k)for(n= 1; n<=4; n++) if(n!=l && n!=k && n!=m){
    		  	 T2[n][m][l][k]=(Q[m][l][k]+P[m][l][k]*F[n][m][l][k]-beta[m][l][k]*E[m][l][k]
					 -sqrt(pow(Q[m][l][k]+P[m][l][k]*F[n][m][l][k]-beta[m][l][k]*E[m][l][k],2.0)-4.0*P[m][l][k]
							*(Q[m][l][k]*F[n][m][l][k] +beta[m][l][k]*ms[k]-I*beta[m][l][k]*rho )) )/(-P[m][l][k]);
			}	


//cout the value of a_{lk} ... F_{nmlk}
	cout << "Khiem's Code" << endl << "q_ij" << endl;
	cout << q10 << "\t" << q11 << "\t" << q12 << "\t" << q13 << "\t" << endl;
	cout << q20 << "\t" << q21 << "\t" << q22 << "\t" << q23 << "\t" << endl;
	cout << q30 << "\t" << q31 << "\t" << q32 << "\t" << q33 << "\t" << endl;
	cout << q40 << "\t" << q41 << "\t" << q42 << "\t" << q43 << "\t" << endl;
	cout << "m_i square" << endl;
	cout << m1s << "\t" << m2s << "\t" << m3s << "\t" << m4s << "\t" << endl;
	
	cout << "terms\treal part\timage part" << endl;
	for(k=1; k<=4; k++){
		for(l=1; l<=4; l++) if(l!=k) {
			cout << "a_" << k << l << "\t" << real_part(a[k][l]) << "\t" << imag_part(a[k][l]) << endl;
		}
	}	
	for(k=1; k<=4; k++) for(l=1; l<=4; l++) if(l!=k) {
		cout << "b_" << k << l << "\t" << real_part(b[k][l]) << "\t" << imag_part(b[k][l]) << endl;
	}		 
	for(k=1; k<=4; k++) for(l=1; l<=4; l++) if(l!=k) {
		cout << "c_" << k << l << "\t" << real_part(c[k][l]) << "\t" << imag_part(c[k][l]) << endl;
	} 
	for(k=1; k<=4; k++) for(l=1; l<=4; l++) if(l!=k) {
		cout << "d_" << k << l << "\t" << real_part(d[k][l]) << "\t" << imag_part(d[k][l]) << endl;
	} 
	for(k=1; k<=4; k++) for(l=1; l<=4; l++) if(l!=k) {
		cout << "AC_" << k << l << "\t" << real_part(AC[k][l]) << "\t" << imag_part(AC[k][l]) << endl;
	} 
	for(k=1; k<=4; k++) for(l=1; l<=4; l++) if(l!=k) {
		cout << "alpha_" << k << l << "\t" << real_part(alpha[k][l]) << "\t" << imag_part(alpha[k][l]) << endl;
	} 
	for(k=1; k<=4; k++) for(l=1; l<=4; l++) if(l!=k) {
		cout << "f_" << k << l << "\t" << real_part(f[k][l]) << "\t" << imag_part(f[k][l]) << endl;
	} 
	for(k=1; k<=4; k++) for(l=1; l<=4; l++) if(l!=k) {
		cout << "fminus_" << k << l << "\t" << real_part(fminus[k][l]) << "\t" << imag_part(fminus[k][l]) << endl;
	} 
	for(k=1; k<=4; k++) for(l=1; l<=4; l++) if(l!=k) for(m=1; m<=4; m++) if(m!=l && m!=k) {
		cout << "A_" << m << k << l << "\t" << real_part(A[m][k][l]) << "\t" << imag_part(A[m][k][l]) << endl;
	} 
	for(k=1; k<=4; k++) for(l=1; l<=4; l++) if(l!=k) for(m=1; m<=4; m++) if(m!=l && m!=k) {
		cout << "B_" << m << k << l << "\t" << real_part(B[m][k][l]) << "\t" << imag_part(B[m][k][l]) << endl;
	} 
	for(k=1; k<=4; k++) for(l=1; l<=4; l++) if(l!=k) for(m=1; m<=4; m++) if(m!=l && m!=k) {
		cout << "C_" << m << k << l << "\t" << real_part(C[m][k][l]) << "\t" << imag_part(C[m][k][l]) << endl;
	} 
	for(k=1; k<=4; k++) for(l=1; l<=4; l++) if(l!=k) for(m=1; m<=4; m++) if(m!=l && m!=k) {
		cout << "D_" << m << k << l << "\t" << real_part(D[m][k][l]) << "\t" << imag_part(D[m][k][l]) << endl;
	} 
	for(k=1; k<=4; k++) for(l=1; l<=4; l++) if(l!=k) for(m=1; m<=4; m++) if(m!=l && m!=k) {
		cout << "beta_" << m << k << l << "\t" << real_part(beta[m][k][l]) << "\t" << imag_part(beta[m][k][l]) << endl;
	} 
	for(k=1; k<=4; k++) for(l=1; l<=4; l++) if(l!=k) for(m=1; m<=4; m++) if(m!=l && m!=k) {
		cout << "phi_" << m << k << l << "\t" << real_part(phi[m][k][l]) << "\t" << imag_part(phi[m][k][l]) << endl;
	} 
	for(k=1; k<=4; k++) for(l=1; l<=4; l++) if(l!=k) for(m=1; m<=4; m++) if(m!=l && m!=k) {
		cout << "Q_" << m << k << l << "\t" << real_part(Q[m][k][l]) << "\t" << imag_part(Q[m][k][l]) << endl;
	}
	for(k=1; k<=4; k++) for(l=1; l<=4; l++) if(l!=k) for(m=1; m<=4; m++) if(m!=l && m!=k) {
		cout << "P_" << m << k << l << "\t" << real_part(P[m][k][l]) << "\t" << imag_part(P[m][k][l]) << endl;
	} 
	for(k=1; k<=4; k++) for(l=1; l<=4; l++) if(l!=k) for(m=1; m<=4; m++) if(m!=l && m!=k) {
		cout << "E_" << m << k << l << "\t" << real_part(E[m][k][l]) << "\t" << imag_part(E[m][k][l]) << endl;
	}
	for(k=1; k<=4; k++) for(l=1; l<=4; l++) for(m=1; m<=4; m++) {
		cout << "g_" << m << k << l << "\t" << real_part(g[m][k][l]) << "\t" << imag_part(g[m][k][l]) << endl;
	} 
	for(k=1; k<=4; k++) for(l=1; l<=4; l++) if(l!=k) for(m=1; m<=4; m++) if(m!=l && m!=k) {
		cout << "gminus_" << m << k << l << "\t" << real_part(gminus[m][k][l]) << "\t" << imag_part(gminus[m][k][l]) << endl;
	} 
	for(k=1; k<=4; k++) for(l=1; l<=4; l++) if (l!=k) for(m=1; m<=4; m++) if(m!=l && m!=k) for(n=1; n<=4; n++) if(n!=l && n!=k && n!=m){
		cout << "F" << n << m << l << k << "\t" << real_part(F[n][m][l][k]) << "\t" << imag_part(F[n][m][l][k]) << endl;
	}
	for(k=1; k<=4; k++) for(l=1; l<=4; l++) if(l!=k) for(m=1; m<=4; m++) if(m!=l && m!=k) {
		cout << "z1beta_" << m << k << l << "\t" << real_part(z1beta[m][k][l]) << "\t" << imag_part(z1beta[m][k][l]) << endl;
	} 
	for(k=1; k<=4; k++) for(l=1; l<=4; l++) if(l!=k) for(m=1; m<=4; m++) if(m!=l && m!=k) {
		cout << "z2beta_" << m << k << l << "\t" << real_part(z2beta[m][k][l]) << "\t" << imag_part(z2beta[m][k][l]) << endl;
	} 
	for(k=1; k<=4; k++) for(l=1; l<=4; l++) if(l!=k) for(m=1; m<=4; m++) if(m!=l && m!=k) {
		cout << "z1phi_" << m << k << l << "\t" << real_part(z1phi[m][k][l]) << "\t" << imag_part(z1phi[m][k][l]) << endl;
	} 
	for(k=1; k<=4; k++) for(l=1; l<=4; l++) if(l!=k) for(m=1; m<=4; m++) if(m!=l && m!=k) {
		cout << "z2phi_" << m << k << l << "\t" << real_part(z2phi[m][k][l]) << "\t" << imag_part(z2phi[m][k][l]) << endl;
	} 
	for(k=1; k<=4; k++) for(l=1; l<=4; l++) if (l!=k) for(m=1; m<=4; m++) if(m!=l && m!=k) for(n=1; n<=4; n++) if(n!=l && n!=k && n!=m){
		cout << "T1" << n << m << l << k << "\t" << real_part(T1[n][m][l][k]) << "\t" << imag_part(T1[n][m][l][k]) << endl;
	}
	for(k=1; k<=4; k++) for(l=1; l<=4; l++) if (l!=k) for(m=1; m<=4; m++) if(m!=l && m!=k) for(n=1; n<=4; n++) if(n!=l && n!=k && n!=m){
		cout << "T2" << n << m << l << k << "\t" << real_part(T2[n][m][l][k]) << "\t" << imag_part(T2[n][m][l][k]) << endl;
	}
	return EXIT_SUCCESS;								
}

// theta function defind as theta(x)= 1 if x>=0; and =0 if x<0;
ex theta(ex x){
	ex temp;
	if(x>=0.0){
		temp=1.0;
	}
	else {
		temp=0.0;
	}
	return temp;
}
// delta function defind as delta(x)=0 if x<>0; and =1 if x==0.
ex delta(ex x){
	ex temp;
	if (x!=0.0){
		temp=0.0;
	}
	else {
		temp=1.0;
	}
	return temp;
}
// define Etafunction here
ex Etafunction (ex x, ex y){
	return (  theta(-imag_part(x)) *theta(-imag_part(y))*theta(imag_part(x*y))
	          - theta(imag_part(x)) *theta(imag_part(y))*theta(-imag_part(x*y))   )*2.0*I*Pi;
}
// define R function here.
ex R(ex x, ex y){
	if(x==y){
		cout<<"R function have denominate zero"<< endl;
	}
	else{
		return  (log(x)-log(y))/(x-y);
	}
}
// define Gpositive here.	
ex Gpositive (ex x,ex y){
	if(x==y){
		cout<<"Gpositive function have denominate zero"<< endl;
	}
	else{
		return R(-x,-y);
	}
}
// define Gnegative here.
ex Gnegative (ex x,ex y){
	if(x==y){
		cout<<"Gnegative function have denominate zero"<< endl;
	}
	else{
		return R(x,y);
	}
}
// define LnGpositive here 
ex LnGpositive(ex x, ex y, ex z){
	ex temp;
	if(x==y){
		cout<<"LnGpositive function have denominate zero"<< endl;
	}
	else {
		temp= R(-x,-y)
			+( log(x)-log(y)-log(x)*log(x)/2.0  +log(y)*log(y)/2.0  + Li2(1.0-z/y)  - Li2(1.0-z/x)  )/(x-y)
			+( log(-y)*( Etafunction(y-z,1.0/(1.0+y)) -Etafunction(y-z,1.0/y))- log(-x)*( Etafunction(x-z,1.0/(1.0+x)) -Etafunction(x-z,1.0/x)))/(x-y)
			+(log(1.0-z/y)*Etafunction(z,1.0/y)-log(1.0-z/x)*Etafunction(z,1.0/x) )/(x-y)  ;		
	}
	return temp;
}
// define LnGnegative here.
ex LnGnegative(ex x, ex y, ex z){
	ex temp;
	if(x==y){
		cout<<"LnGnegative function have denominate zero"<< endl;
	}
	else {
		temp= (-log(z)+log(-z)+1.0)*R(x,y) 
			-( log(x)-log(y)-log(x)*log(x)/2.0+log(y)*log(y)/2.0 + Li2(1.0-z/y)-Li2(1.0-z/x) )/(x-y)	
			-(log(y)*( Etafunction(z-y,1.0/(1.0-y))-Etafunction(z-y,-1.0/y) )-log(x)*(Etafunction(z-x,1.0/(1.0-x)) -Etafunction(z-x,-1.0/x) ) )/(x-y)
			- ( log(1.0-z/y)*Etafunction(z,1.0/y)-log(1.0-z/x)*Etafunction(z,1.0/x) )/(x-y) ;
	}
	return temp;
	
}


 void input_option(int argc, char *argv[], lst &p, lst &m){
	//Use optget_long to get the long options with their values
	ex m1, m2, m3, m4;
	m1 = symbol("m1"); m2 = symbol("m2"); m3 = symbol("m3"); m4 = symbol("m4");
	
	// Store the P's for loop tool
	ex p1s=symbol("p1s"),p2s=symbol("p2s"),p3s=symbol("p3s"),p4s=symbol("p4s"),p12s=symbol("p12s"),p23s=symbol("p23s");
	
	
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
			{"p1s", 1, 0, 0},//P1 Square
			{"p2s", 1, 0, 0},//P2 Square
			{"p3s", 1, 0, 0},//P3 Square
			{"p4s", 1, 0, 0},//P4 Square
			{"p12s", 1, 0, 0},//(P1+P2) Square
			{"p23s", 1, 0, 0},//(P2+P3) Square
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
			case 4: p1s = temp_numcast; 
				loopToolInput_Indicator = 1;
				break;
			case 5: p2s = temp_numcast; 
				loopToolInput_Indicator = 1;
				break;
			case 6: p3s = temp_numcast; 
				loopToolInput_Indicator = 1;
				break;
			case 7: p4s = temp_numcast; 
				loopToolInput_Indicator = 1;
				break;
			case 8: p12s = temp_numcast; 
				loopToolInput_Indicator = 1;
				break;
			case 9: p23s = temp_numcast; 
				loopToolInput_Indicator = 1;
				break;
			}
			break;
		default:
			throw std::runtime_error("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nUnknow input option, please use:\n --m1s --m2s --m3s --m4s --q10 --q20 --q21 --q30 --q31 --q32\nOR\n--m1s --m2s --m3s --m4s --p1s --p2s --p3s --p4s --p12s --p23s\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
		}
	}
		//input m2s, m3s, m4s, m1s for xloop
	m = lst(m2, m3, m4, m1);
	p = lst(p1s, p2s, p3s, p4s, p12s, p23s);
}

// end.........the end......
