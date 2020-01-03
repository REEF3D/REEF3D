/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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
#include"fillvec.h"
#include"weno_nug_func.h"

class flux;
class cpt;

#ifndef WENO_FLUX_N_NUG_H_
#define WENO_FLUX_N_NUG_H_

using namespace std;

class weno_flux_N_nug : public convection, public weno_nug_func, public fillvec
{
public:
	weno_flux_N_nug(lexer*);
	virtual ~weno_flux_N_nug();

	virtual void start(lexer*,fdm*,field&,int,field&,field&,field&);

private:
    double aij(lexer*, fdm*, field&, int,field&,field&,field&,double*,double*,double*, cpt&);
    
    double aij_sig(lexer*, fdm*, field&, int,field&,field&,field&,double*,double*,double*, cpt&);

	virtual double fx(lexer*, fdm*, field&, field&, int, double, cpt&);
	virtual double fy(lexer*, fdm*, field&, field&, int, double, cpt&);
	virtual double fz(lexer*, fdm*, field&, field&, int, double, cpt&);
	void iqmin(lexer*, fdm*, field&, field&, cpt&);
	void jqmin(lexer*, fdm*, field&, field&, cpt&);
	void kqmin(lexer*, fdm*, field&, field&, cpt&);
	void iqmax(lexer*, fdm*, field&, field&, cpt&);
	void jqmax(lexer*, fdm*, field&, field&, cpt&);
	void kqmax(lexer*, fdm*, field&, field&, cpt&);


	double L,grad;

	double gradx, grady, gradz;
	double fu1,fv1,fw1,fu2,fv2,fw2;
    
    double ivel1,ivel2,jvel1,jvel2,kvel1,kvel2;

    
    flux *pflux;
    

};

#endif
