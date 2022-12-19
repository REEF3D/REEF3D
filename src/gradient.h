/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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

class fdm;
class field;
class lexer;
class slice;

#ifndef GRADIENT_H_
#define GRADIENT_H_

using namespace std;

class gradient : virtual public increment
{
public:

	gradient(lexer*);
	 ~gradient();



	//x
	 double xdx(fdm*, field&);
	 double xdy(fdm*, field&);
	 double xdz(fdm*, field&);

	 double xdxx(fdm*, field&);
	 double xdyy(fdm*, field&);
	 double xdzz(fdm*, field&);

	 double xdxy(fdm*, field&);
	 double xdxz(fdm*, field&);
	 double xdyz(fdm*, field&);
	//y

	 double ydx(fdm*, field&);
	 double ydy(fdm*, field&);
	 double ydz(fdm*, field&);

	 double ydxx(fdm*, field&);
	 double ydyy(fdm*, field&);
	 double ydzz(fdm*, field&);

	 double ydxy(fdm*, field&);
	 double ydxz(fdm*, field&);
	 double ydyz(fdm*, field&);

	//z
	 double zdx(fdm*, field&);
	 double zdy(fdm*, field&);
	 double zdz(fdm*, field&);

	 double zdxx(fdm*, field&);
	 double zdyy(fdm*, field&);
	 double zdzz(fdm*, field&);

	 double zdxy(fdm*, field&);
	 double zdxz(fdm*, field&);
	 double zdyz(fdm*, field&);

	 double ddwenox(fdm*, field&, double);
	 double ddwenoy(fdm*, field&, double);
	 double ddwenoz(fdm*, field&, double);
	 
	 double dfwenox(fdm*, field&, double);
	 double dfwenoy(fdm*, field&, double);
	 double dfwenoz(fdm*, field&, double);
     
     double dslwenox(fdm*, slice&, double);
	 double dslwenoy(fdm*, slice&, double);

	//--------------------------------

	//p
	 double pudx(lexer*,fdm*);
	 double pudy(lexer*,fdm*);
	 double pudz(lexer*,fdm*);

	 double pvdx(lexer*,fdm*);
	 double pvdy(lexer*,fdm*);
	 double pvdz(lexer*,fdm*);

	 double pwdx(lexer*,fdm*);
	 double pwdy(lexer*,fdm*);
	 double pwdz(lexer*,fdm*);
     
     //
     double pudx(lexer*,field&);
	 double pudy(lexer*,field&);
	 double pudz(lexer*,field&);

	 double pvdx(lexer*,field&);
	 double pvdy(lexer*,field&);
	 double pvdz(lexer*,field&);

	 double pwdx(lexer*,field&);
	 double pwdy(lexer*,field&);
	 double pwdz(lexer*,field&);
	 

	//--------------------------------

	//u
	 double udx(fdm*);
	 double udy(fdm*);
	 double udz(fdm*);

	 double udxx(fdm*);
	 double udyy(fdm*);
	 double udzz(fdm*);

	//v
	 double vdx(fdm*);
	 double vdy(fdm*);
	 double vdz(fdm*);

	 double vdxx(fdm*);
	 double vdyy(fdm*);
	 double vdzz(fdm*);

	//w
	 double wdx(fdm*);
	 double wdy(fdm*);
	 double wdz(fdm*);

	 double wdxx(fdm*);
	 double wdyy(fdm*);
	 double wdzz(fdm*);
	 

	double grad1,grad2;
	double grad;
	const double dx;

	void iqmin(field&, double);
	void jqmin(field&, double);
	void kqmin(field&, double);
	void iqmax(field&, double);
	void jqmax(field&, double);
	void kqmax(field&, double);
    
    void iqminsl(slice&, double);
	void jqminsl(slice&, double);
	void iqmaxsl(slice&, double);
	void jqmaxsl(slice&, double);

	const double tttw,fourth,third,sevsix,elvsix,sixth,fivsix,tenth;
	const double sixten,treten;
	const double epsilon,smallnum;
	double is1,is2,is3;
	double alpha1,alpha2,alpha3;
	double w1,w2,w3;
	double q1,q2,q3,q4,q5;
	double gradx, grady, gradz;
	double f1,f2,f3,f4;
	
	double ddx,ddy,ddz;
	double fu1,fv1,fw1,fu2,fv2,fw2;
	double fx(fdm*, field&, double);
	double fy(fdm*, field&, double);
	double fz(fdm*, field&, double);
	void iqfmin(field&);
	void jqfmin(field&);
	void kqfmin(field&);
	void iqfmax(field&);
	void jqfmax(field&);
	void kqfmax(field&);
	void isf(field&);
	void alphaf();
	void weightf();

	void is();
	void alpha();
	void weight();
    
    lexer *p;
};

#endif
