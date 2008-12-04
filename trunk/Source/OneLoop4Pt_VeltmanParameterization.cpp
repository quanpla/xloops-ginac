/*******************************************************************************
**
**	XLOOP-GiNaCs Project
**		1Loop4Pt Veltman Parameterization using Vegas Calculus Integral.
**		
**
**	HCMUNS,
**	November 2008
**
**
**	Author(s): 	Son, Do (sondo01@gmail.com)
**				Khiem, Phan (phanhongkhiem@gmail.com)
**				Quan, Phan (anhquan.phanle@gmail.com)
********************************************************************************
**
**	Notes:
**
********************************************************************************
**
** Historial Log:
**	Date		Version	Author		Description
**	_________	_______	_________	________________________________
**	20081102	1.0		Quan Phan	Realign and create the declaration file,
**									using GiNaC to calculate
**	200810..	1.0		Khiem Phan	Create this file
*******************************************************************************/

/*
 *  xloops Copyright (C) 1997,2000-2002 Johannes Gutenberg University Mainz,
 *  Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */


#include <iostream>
#include <sstream>

#include <getopt.h>


extern "C" {
#include "vegas.h" 
}

#include "OneLoop4Pt_VeltmanParameterization.h"

using namespace std;
using namespace GiNaC;
using namespace XLOOPS_1L4Pt_VeltmanParameterization;

int main(int argc, char **argv){
	ex m, p;
	getInput(argc, argv);
	preCalculate();
	
	double real_result, real_error, imag_result, imag_error;
	
	callVegas(1, real_result, real_error); // 1 = real case, 0 = image case
	callVegas(0, imag_result, imag_error);
	
	cout << "\n*****************************************************************************\n"
		<<   "**  (imag, real) D0 =    (" << real_result << " +/- " << real_error << "  ,    " << imag_result << " +/- " << imag_error
		<<   "*****************************************************************************\n";
	return EXIT_SUCCESS;
}



namespace XLOOPS_1L4Pt_VeltmanParameterization{

void preCalculate(){
/*
	Precalculate D0 to be substituted by running variables in the integrand
*/
	//	1.	Read the input p & m
	ex p1s = p.op(0), p2s = p.op(1), p3s = p.op(2), p4s = p.op(3), p12s = p.op(4), p23s = p.op(5);
	ex m1s = m.op(0), m2s = m.op(1), m3s = m.op(2), m4s = m.op(3);
	
	//	2.	Calculate intermidate terms
	ex 	f1 = -(p12s+p1s-p2s), f2 = -(-p23s+p1s+p4s), f3 = -(p12s+p2s-p3s),
		f4 =  (p1s-m2s+m1s),  f5 =  (p12s-m3s+m1s),  f6 =  (p4s-m4s+m1s);
	
	//	3.	Calculate it
	D0 = -p1s * pow(x0,2)
		-p12s * pow(1-x0,2) * pow(x1,2)
		-p4s * pow(1-x0,2) * pow(1-x1,2) * pow(x2,2) 
		+f1 * x0 * (1-x0) * x1
		+f2 * x0 * (1-x0) * x2 * (1-x1)
		+f3 * pow(1-x0,2) * (1-x1) * x1 * x2
		+f4 * x0
		+f5 * (1-x0) * x1
		+f6 * (1-x0) * (1-x1) * x2
		-m1s;
	
	D0 = pow(1-x0,2) * (1-x1) / pow(D0, 2);
	//evalf it to make it faster substitute;
	D0 = D0.evalf();
}

void pr_val(double x[DIMENSION], double f[FUNCTIONS])
/*
	This is the integrand
*/
{
	//	1.	Substitute the values into the D0 expression
	ex cur_D0 = D0.subs(lst(x0 == x[0], x1 == x[1], x2 == x[2]));
	cur_D0 = cur_D0.evalf();
	if (is_a<numeric>(cur_D0)){
		// check if the calculate is a numeric, else raise error
		if (is_RealCalculation == 1)
			f[0] = ex_to<numeric>(real_part(cur_D0)).to_double();
		else
			f[0] = ex_to<numeric>(imag_part(cur_D0)).to_double();
	}
	else{
		stringstream err_msg;
		err_msg << "Cannot evaluation D0 to float. D0 = " << cur_D0 << endl;
		throw runtime_error(err_msg.str());
	}
}

void callVegas(int is_RealCalculation_, double &value, double &error)
{
	//	Declarations
	int i;
	
	double estim[FUNCTIONS];   /* estimators for integrals                     */
	double std_dev[FUNCTIONS]; /* standard deviations                          */
	double chi2a[FUNCTIONS];   /* chi^2/n                                      */
	double reg[2*DIMENSION];   /* integration domain                           */
	
	int vegas_output_option = NPRN_INPUT | NPRN_RESULT;
	//	1.	Initialize
	for (i = 0; i < DIMENSION; i++) {
		reg[i] = 0.0;
		reg[i+DIMENSION] = 1.0;
	}
	
	// decide to calculate real part or image part
	is_RealCalculation = is_RealCalculation_;
	
	vegas(reg, DIMENSION, pr_val, 0, numOfSample_s1, numOfIterations, vegas_output_option, FUNCTIONS, 0, 1, estim, std_dev, chi2a);
	vegas(reg, DIMENSION, pr_val, 1, numOfSample_s2, numOfIterations, vegas_output_option, FUNCTIONS, 0, 2, estim, std_dev, chi2a);
	
	// 
	value = estim[0];
	error = std_dev[0]; 
}

void getInput(int argc, char *argv[]){
	/*Use optget_long to get the long options with their values*/
	ex m1, m2, m3, m4,
		p1s, p2s, p3s, p4s, p12s, p23s;
	
	m1 = symbol("m1s"); m2 = symbol("m2s"); m3 = symbol("m3s"); m4 = symbol("m4s");
	
	p1s=symbol("p1s"); p2s=symbol("p2s"); p3s=symbol("p3s"); p4s=symbol("p4s");
	p12s=symbol("p12s"); p23s=symbol("p23s");
	
	int c;
	int digit_optind = 0;
	double temp_numcast = 0;
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
				break;
			case 5: p2s = temp_numcast; 
				break;
			case 6: p3s = temp_numcast; 
				break;
			case 7: p4s = temp_numcast; 
				break;
			case 8: p12s = temp_numcast; 
				break;
			case 9: p23s = temp_numcast; 
				break;
			}
			break;
		default:
			throw std::runtime_error("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nUnknow input option, please use:\n--m1s --m2s --m3s --m4s --p1s --p2s --p3s --p4s --p12s --p23s\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
		}
	}
	// there is a rotation b/c of Khiem's transformation
	m = lst(m4, m1, m2, m3);
	p = lst(p1s, p2s, p3s, p4s, p12s, p23s);
}

} // namespace XLOOPS_1L4Pt_VeltmanParameterization
