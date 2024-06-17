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

#ifndef NHFLOW_SCALAR_IWENO_H_
#define NHFLOW_SCALAR_IWENO_H_

#include"nhflow_scalar_convection.h"
#include"weno_nug_func.h"

class nhflow_scalar_advec;

using namespace std;

class nhflow_scalar_iweno : public nhflow_scalar_convection, public weno_nug_func
{
public:
	nhflow_scalar_iweno (lexer*);
	virtual ~nhflow_scalar_iweno();

	virtual void start(lexer*,fdm_nhf*,double*,int,double*,double*,double*);

private:

	void aij(fdm_nhf*, double*,double*,double*,double*);
	void aij_south(lexer*,fdm_nhf*,double*, double*);
	void aij_north(lexer*,fdm_nhf*,double*, double*);
	void aij_east(lexer*,fdm_nhf*,double*, double*);
	void aij_west(lexer*,fdm_nhf*,double*, double*);
	void aij_top(lexer*,fdm_nhf*,double*, double*);
	void aij_bottom(lexer*,fdm_nhf*,double*, double*);
    
    void iqmin(lexer*, fdm_nhf*, double*);
	void jqmin(lexer*, fdm_nhf*, double*);
	void kqmin(lexer*, fdm_nhf*, double*);
	void iqmax(lexer*, fdm_nhf*, double*);
	void jqmax(lexer*, fdm_nhf*, double*);
	void kqmax(lexer*, fdm_nhf*, double*);

	const double tttw,fourth,third,sevsix,elvsix,sixth,fivsix,tenth;
	const double sixten,treten,epsi;


	double is1,is2,is3;
	double alpha1,alpha2,alpha3;
	double w1,w2,w3;
	double umin, umax, uplus;
	int count;

    
    
    double ivel1,ivel2,jvel1,jvel2,kvel1,kvel2;
    double iadvec,jadvec,kadvec;
    
    nhflow_scalar_advec *padvec;
    
    double *DX,*DY,*DZ;
    
};

#endif
