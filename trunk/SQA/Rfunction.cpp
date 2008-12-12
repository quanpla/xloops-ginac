#include <iostream>
#include "ginac/ginac.h"
#include "ginac/flags.h"
#include <sstream>
extern "C" {
#include "vegas.h"
}
#include <stdexcept>
#include <fstream>
using namespace std;
using namespace GiNaC;

#define DIMENSION 1
#define FUNCTIONS 1
ex theta(ex );
ex delta(ex );
ex Etafunction(ex ,ex );
static int imagCompute_flag; // =0: compute real part, = 1: compute imaginary part
//for(double zz = startZ; zz <= endZ; zz += stepZ){
ex x1,y;
ex a,b,c,d;

void val(double x[DIMENSION], double f1[FUNCTIONS])
{
	
	ex D0;
	D0 = (Pi/2.0)*(1.0+tan(Pi*x[0]/2.0)*tan(Pi*x[0]/2.0))/ ( (tan(Pi*x[0]/2.0)+x1)*(tan(Pi*x[0]/2.0)+y) );
	ex out_D0=D0.evalf();
	//separate real and image part
        if (imagCompute_flag == 0){
		out_D0 = real_part(out_D0).evalf();
		double num = ex_to<numeric>(out_D0).to_double();
		f1[0] = num;
	}
	else{
		out_D0 = imag_part(out_D0).evalf();
		double num = ex_to<numeric>(out_D0).to_double();
		f1[0] = num;
	}	
 } 


int main(int argc, char *argv[]){
	for( b= -50.0; b <= 50.0; b += 1.0)
	for( c= -50.0; c <= 50.0; c += 1.0) if( b!=c ){
			x1=10.0+I*b;
			y=c+I*10.0;
	double estim[FUNCTIONS]; /* estimators for integrals */
	double std_dev[FUNCTIONS]; /* standard deviations */
	double chi2a[FUNCTIONS]; /* chi^2/n */
	double reg[2*DIMENSION]; /* integration domain */
	//	1.	Initializing
	for (int i=0; i<DIMENSION; i++) {
		reg[i] = 0.0;
		reg[i+DIMENSION] =1.0;
	}

	int 	number_of_iterations =10,
		number_of_loops_per_iteration1 = 1000,
		number_of_loops_per_iteration2 = 10000;

	//cout << endl << "*******************************************************************************" << endl;
	//cout << "Starting computing Real part......" << endl;
	imagCompute_flag = 0; // compute the real part
	
	// 1st stage
	//cout << endl << "Real part 1st stage....." << endl;
	vegas(reg, DIMENSION, val, 0, number_of_loops_per_iteration1, number_of_iterations, 0/*No intermidate output NPRN_INPUT | NPRN_RESULT*/, 1, 0, 1, estim, std_dev, chi2a);
	
	// 2nd stage
	//cout << endl << "Real part 2nd stage....." << endl;
	vegas(reg, DIMENSION, val, 1, number_of_loops_per_iteration2, number_of_iterations, 0/*No intermidate output NPRN_INPUT | NPRN_RESULT*/, FUNCTIONS, 0, 2, estim, std_dev, 	chi2a);
	
	// get the result
	double real_part = estim[0];
	double error = std_dev[0];
	
	
	
	//cout << endl << endl << "*******************************************************************************" << endl;
	//cout << "Starting computing Imaginary part......" << endl;
	imagCompute_flag = 1 ; //Compute the Imaginary part
	// 1st stage
	//cout << endl << "Imaginary part 1st stage....." << endl;
	vegas(reg, DIMENSION, val, 0, number_of_loops_per_iteration1, number_of_iterations, 0/*No intermidate output NPRN_INPUT | NPRN_RESULT*/, 1, 0, 1, estim, std_dev, chi2a);
	
	// 2nd stage
	//cout << endl << "Imaginary part 2nd stage....." << endl;
	vegas(reg, DIMENSION, val, 1, number_of_loops_per_iteration2, number_of_iterations, 0/*No intermidate output NPRN_INPUT | NPRN_RESULT*/, FUNCTIONS, 0, 2, estim, std_dev, chi2a);
	
	// get the result
	double im_part = estim[0];
	double error2 = std_dev[0]; 
	
	// QP reformat to show only necessary info.
	//cout << endl<< endl << "*************************************************************************************" << endl;
	//cout << "**  Final result (real +/- error, imaginary +/- error):" << endl << "**\t\t";
	cout << real_part << "\t" << error << "\t" << im_part << "\t" << error2 <<endl;
	// the ginac result
	
	ex GiNaCRfunction= (log(x1)-log(y))/(x1-y);
	cout << real_part(x1) << "\t" << imag_part(x1) << "\t" << real_part(y) << "\t" << imag_part(y) << "\t" << real_part << "\t" << error << "\t" << im_part << "\t" << error2 << "\t" << GiNaCRfunction.evalf() << endl;
	cout << (GiNaCRfunction).evalf() << endl;
	cout << x1<<"and "<< y<<endl;
	}
 return EXIT_SUCCESS;								
}
// delta function defind as delta(x)=0 if x<>0; and =1 if x==0.
ex theta(ex a)
{
if (a>=0.0)
   { 
   return 1.0; 
   }
else 
  {
   return 0.0;
   }
}
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
ex Etafunction(ex x,ex y){
	ex temp;
	temp=2.0*Pi*I*( theta(-imag_part(x))*theta(-imag_part(y))*theta(imag_part(x*y))
	          -theta(imag_part(x))*theta(imag_part(y))*theta(-imag_part(x*y))
                     ) ;
	return temp;
}
