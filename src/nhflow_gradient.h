/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#ifndef NHFLOW_GRADIENT_H_
#define NHFLOW_GRADIENT_H_

using namespace std;

class nhflow_gradient : virtual public increment
{
public:

	nhflow_gradient(lexer*);
	 ~_nhflow_gradient();


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
};

#endif
