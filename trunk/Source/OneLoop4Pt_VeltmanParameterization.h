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
**	20081102	1.0	Quan Phan	Realign and create the declaration file,
**							using GiNaC to calculate
**	200810..	1.0	Khiem Phan	Create this file
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

#ifndef __XLOOPS_ONELOOP_4PT_VeltmanParam_H__
#define __XLOOPS_ONELOOP_4PT_VeltmanParam_H__


// Include GiNaC for precise algebra
#include <ginac/ginac.h>

// Define the size of the matrix of Variables and Functions to put into VEGAS calculations
#define DIMENSION 3
#define FUNCTIONS 1

using namespace GiNaC;

namespace XLOOPS_1L4Pt_VeltmanParameterization{
	// Vegas related Variables and Functions
int 	numOfIterations = 100, // number of iterations in Vegas
	numOfSample_s1 = 1000, // number of seedings in step 1
	numOfSample_s2 = 10000; //number of seedings in step 2

static int is_RealCalculation = 1; // =0 if calculate Image Part, =1 if calculate Real Part
void pr_val(double x[DIMENSION], double f[FUNCTIONS]); // this procedure return the value of the integrand

	// private variables and functions
ex D0; // D0 expression, precalculate, and later will be substituted
ex m, p; // inputs p & m

	// substitute variables
ex 	x0 = symbol("x0"), x1 = symbol("x1"), x2 = symbol("x2"); // running variables

	// the procedures to grab the inputs
void getInput(int argc, char *argv[]);
	// prepare the calculations that will be substitute in the integrand
void preCalculate();
	// call the Vegas to calculate what we need, return value and error
void callVegas(int is_RealCalculation_, double &value, double &error);

}// namespace XLOOPS_1L4Pt_VeltmanParameterization

#endif // __XLOOPS_ONELOOP_4PT_VeltmanParam_H__
