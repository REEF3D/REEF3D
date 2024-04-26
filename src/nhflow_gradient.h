/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"increment.h"

class fdm_nhf;
class lexer;
class slice;

#ifndef NHFLOW_GRADIENT_H_
#define NHFLOW_GRADIENT_H_

using namespace std;

class nhflow_gradient : virtual public increment
{
public:

	nhflow_gradient(lexer*);
	 ~nhflow_gradient();

    double sx(slice&);
    double sy(slice&);
    double limiter(double, double);
    
    double sxx(slice&);
    double syy(slice&);
    
    double dslwenox(slice&, double);
    double dslwenoy(slice&, double);
    
    void iqminsl(slice&, double);
	void jqminsl(slice&, double);
	void iqmaxsl(slice&, double);
	void jqmaxsl(slice&, double);
	//--------------------------------
    
    double dwenox(double*, double);
    double dwenoy(double*, double);
    
    void iqmin(double*, double);
	void jqmin(double*, double);
	void iqmax(double*, double);
	void jqmax(double*, double);
    //--------------------------------

	//u
	 double dudx(double*);
	 double dudy(double*);
	 double dudz(double*);

	 double dudxx(double*);
	 double dudyy(double*);
	 double dudzz(double*);

	//v
	 double dvdx(double*);
	 double dvdy(double*);
	 double dvdz(double*);

	 double dvdxx(double*);
	 double dvdyy(double*);
	 double dvdzz(double*);

	//w
	 double dwdx(double*);
	 double dwdy(double*);
	 double dwdz(double*);

	 double dwdxx(double*);
	 double dwdyy(double*);
	 double dwdzz(double*);
	 

	double grad1,grad2;
	double grad;
	
    
    lexer *p;
    
private:
    double dfdx_min, dfdx_plus, dfdy_min, dfdy_plus, dfdz_min, dfdz_plus;
    double denom,val;
    
    const double tttw,fourth,third,sevsix,elvsix,sixth,fivsix,tenth;
	const double sixten,treten;
	const double epsilon,smallnum,dx;
	double is1,is2,is3;
	double alpha1,alpha2,alpha3;
	double w1,w2,w3;
	double q1,q2,q3,q4,q5;
	double gradx, grady, gradz;
	double f1,f2,f3,f4;
    
    void is();
	void alpha();
	void weight();
};

#endif
