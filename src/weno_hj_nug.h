/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
#include"weno_nug_func.h"

class flux;

#ifndef WENO_HJ_NUG_H_
#define WENO_HJ_NUG_H_

using namespace std;

class weno_hj_nug : public convection, public weno_nug_func
{
public:
	weno_hj_nug(lexer*);
	virtual ~weno_hj_nug();

	virtual void start(lexer*,fdm*,field&,int,field&,field&,field&);

private:
    double aij(lexer*, fdm*, field&, int,field&,field&,field&,double*,double*,double*);
    
	virtual double fx(lexer*, fdm*, field&, field&, int, double);
	virtual double fy(lexer*, fdm*, field&, field&, int, double);
	virtual double fz(lexer*, fdm*, field&, field&, int, double);
	void iqmin(lexer*, fdm*, field&, field&, int);
	void jqmin(lexer*, fdm*, field&, field&, int);
	void kqmin(lexer*, fdm*, field&, field&, int);
	void iqmax(lexer*, fdm*, field&, field&, int);
	void jqmax(lexer*, fdm*, field&, field&, int);
	void kqmax(lexer*, fdm*, field&, field&, int);


	double L,grad;

	double gradx, grady, gradz;
	double fu1,fv1,fw1,fu2,fv2,fw2;
    
    double ivel1,ivel2,jvel1,jvel2,kvel1,kvel2;
    double iadvec,jadvec,kadvec;
    
    double *DX,*DY,*DZ;
    
    flux *pflux;
    
    //---
    void is();
	void alpha();
	void weight();
    
    double is1,is2,is3;
	double alpha1,alpha2,alpha3;
	double w1,w2,w3;
    
    const double tttw,fourth,third,sevsix,elvsix,sixth,fivsix,tenth;
	const double sixten,treten;
	const double epsilon,smallnum;
};

#endif
