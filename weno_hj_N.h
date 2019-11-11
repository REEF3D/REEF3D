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
class vec;

#ifndef WENO_HJ_N_H_
#define WENO_HJ_N_H_

using namespace std;

class weno_hj_N : public fillvec
{
public:
	weno_hj_N(lexer*);
	virtual ~weno_hj_N();

	virtual void start(lexer*,fdm*,fieldint&,vec&,int,field&,field&,field&);

private:
    double aij(lexer*, fdm*, vec&, fieldint&, int,field&,field&,field&, cpt&);
    
    double aij_sig(lexer*, fdm*, field&, int,field&,field&,field&,double*,double*,double*);

	virtual double ddx(lexer*, fdm*, vec&, cpt&);
	virtual double ddy(lexer*, fdm*, vec&, cpt&);
	virtual double ddz(lexer*, fdm*, vec&, cpt&);
	void iqmin(fdm*,vec&, double, cpt&);
	void jqmin(fdm*,vec&, double, cpt&);
	void kqmin(fdm*,vec&, double, cpt&);
	void iqmax(fdm*,vec&, double, cpt&);
	void jqmax(fdm*,vec&, double, cpt&);
	void kqmax(fdm*,vec&, double, cpt&);

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
    int *sizeM;
    
    flux *pflux;
};

#endif
