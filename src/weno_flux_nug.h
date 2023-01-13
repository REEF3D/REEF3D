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

#ifndef WENO_FLUX_NUG_H_
#define WENO_FLUX_NUG_H_

using namespace std;

class weno_flux_nug : public convection, public weno_nug_func
{
public:
	weno_flux_nug(lexer*);
	virtual ~weno_flux_nug();

	virtual void start(lexer*,fdm*,field&,int,field&,field&,field&);

private:
    double aij(lexer*, fdm*, field&, int,field&,field&,field&,double*,double*,double*);
    
	virtual double fx(lexer*, fdm*, field&, field&, int, double);
	virtual double fy(lexer*, fdm*, field&, field&, int, double);
	virtual double fz(lexer*, fdm*, field&, field&, int, double);
	void iqmin(lexer*, field&, field&, int);
	void jqmin(lexer*, field&, field&, int);
	void kqmin(lexer*, field&, field&, int);
	void iqmax(lexer*, field&, field&, int);
	void jqmax(lexer*, field&, field&, int);
	void kqmax(lexer*, field&, field&, int);

    
	double L,grad;
    double sig;

	double gradx, grady, gradz;
	double fu1,fv1,fw1,fu2,fv2,fw2;
    
    double ivel1,ivel2,jvel1,jvel2,kvel1,kvel2;

    
    flux *pflux;

};

#endif
