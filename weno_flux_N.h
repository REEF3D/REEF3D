/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"convection.h"
#include"fillvec.h"

class flux;
class cpt;

#ifndef WENO_FLUX_N_H_
#define WENO_FLUX_N_H_

using namespace std;

class weno_flux_N : public convection, public fillvec
{
public:
	weno_flux_N(lexer*);
	virtual ~weno_flux_N();

	virtual void start(lexer*,fdm*,field&,int,field&,field&,field&);

private:
    double aij(lexer*, fdm*, field&, int,field&,field&,field&, cpt&);
    
    double aij_sig(lexer*, fdm*, field&, int,field&,field&,field&,double*,double*,double*,cpt&);

	virtual double fx(lexer*, fdm*, field&, field&, double, cpt&);
	virtual double fy(lexer*, fdm*, field&, field&, double, cpt&);
	virtual double fz(lexer*, fdm*, field&, field&, double, cpt&);
	void iqmin(lexer*, fdm*, field&, field&, cpt&);
	void jqmin(lexer*, fdm*, field&, field&, cpt&);
	void kqmin(lexer*, fdm*, field&, field&, cpt&);
	void iqmax(lexer*, fdm*, field&, field&, cpt&);
	void jqmax(lexer*, fdm*, field&, field&, cpt&);
	void kqmax(lexer*, fdm*, field&, field&, cpt&);

	double L,grad;
	const double tttw,fourth,third,sevsix,elvsix,sixth,fivsix,tenth;
	const double sixten,treten;
	const double epsilon,smallnum;
	double is1,is2,is3;
	double alpha1,alpha2,alpha3;
	double w1,w2,w3;
	double q1,q2,q3,q4,q5;
	double gradx, grady, gradz;
	double fu1,fv1,fw1,fu2,fv2,fw2;
    int n0;
    
    double ivel1,ivel2,jvel1,jvel2,kvel1,kvel2;


	void is(field&);
	void alpha();
	void weight();
    
    flux *pflux;
};

#endif
