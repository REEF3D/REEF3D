/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"convection.h"
#include"increment.h"

class flux;
class cpt;

#ifndef WENO_HJ_H_
#define WENO_HJ_H_

using namespace std;

class weno_hj : public convection, public increment
{
public:
	weno_hj(lexer*);
	virtual ~weno_hj();

	virtual void start(lexer*,fdm*,field&,int,field&,field&,field&);

private:
    double aij(lexer*, fdm*, field&, int,field&,field&,field&,double*,double*,double*);
    
	virtual double ddx(lexer*, fdm*, field&);
	virtual double ddy(lexer*, fdm*, field&);
	virtual double ddz(lexer*, fdm*, field&);
	void iqmin(field&, double, double*);
	void jqmin(field&, double, double*);
	void kqmin(field&, double, double*);
	void iqmax(field&, double, double*);
	void jqmax(field&, double, double*);
	void kqmax(field&, double, double*);

	double L,grad;
	const double tttw,fourth,third,sevsix,elvsix,sixth,fivsix,tenth;
	const double sixten,treten;
	const double epsilon,smallnum;
	double is1,is2,is3;
	double alpha1,alpha2,alpha3;
	double w1,w2,w3;
	double q1,q2,q3,q4,q5;
	double gradx, grady, gradz;
    
    double ivel1,ivel2,jvel1,jvel2,kvel1,kvel2;
    double iadvec,jadvec,kadvec;


	void is();
	void alpha();
	void weight();
    
    flux *pflux;
};

#endif
