/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs
 * 
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

#ifndef IWENO_FLUX_H_
#define IWENO_FLUX_H_

using namespace std;

class iweno_flux : public convection, public increment
{
public:
	iweno_flux (lexer*);
	virtual ~iweno_flux();

	virtual void start(lexer*,fdm*,field&,int,field&,field&,field&);

private:
    void wenoloop1(lexer*,fdm*,field&,int,field&,field&,field&);
    void wenoloop2(lexer*,fdm*,field&,int,field&,field&,field&);
    void wenoloop3(lexer*,fdm*,field&,int,field&,field&,field&);
    void wenoloop4(lexer*,fdm*,field&,int,field&,field&,field&);
	
	void iqmin(fdm*, field&, field&, int);
	void jqmin(fdm*, field&, field&, int);
	void kqmin(fdm*, field&, field&, int);
	void iqmax(fdm*, field&, field&, int);
	void jqmax(fdm*, field&, field&, int);
	void kqmax(fdm*, field&, field&, int);

    void is_1(field&);
	void is_2(field&);

	void alpha_calc_1();
	void alpha_calc_2();

	void weight_calc_1();
	void weight_calc_2();

	void aij(fdm*, field&,field&,field&,field&);
	void aij_x(lexer*,fdm*,field&,field&);
	void aij_y(lexer*,fdm*,field&,field&);
	void aij_z(lexer*,fdm*,field&,field&);

	const double tttw,fourth,third,sevsix,elvsix,sixth,fivsix,tenth;
	const double sixten,treten,epsilon,deltin;


	double is1_1,is2_1,is3_1,is1_2,is2_2,is3_2;
	double alpha1_1,alpha2_1,alpha3_1,alpha1_2,alpha2_2,alpha3_2;
	double w1_1,w2_1,w3_1,w1_2,w2_2,w3_2;
	double umin, umax, uplus;
	double a1,a2,b1,b2,c1,c2;
	double q0,q1,q2,q3,q4,q5,q6;
	int count,rocount,countN,coliN,aiipos;
	int *range;
    
    double ivel1,ivel2,jvel1,jvel2,kvel1,kvel2;
    
    flux *pflux;
};

#endif
