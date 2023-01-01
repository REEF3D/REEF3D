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

#include"nhflow_convection.h"
#include"weno_nug_func.h"

class nhflow_flux;

#ifndef NHFLOW_WENO_FLUX_H_
#define NHFLOW_WENO_FLUX_H_

using namespace std;

class nhflow_weno_flux : public nhflow_convection, public weno_nug_func
{
public:
	nhflow_weno_flux(lexer*);
	virtual ~nhflow_weno_flux();

	virtual void start(lexer*, fdm_nhf*, double*, int, double*, double*,double*);

private:
    double aij(lexer*, fdm_nhf*, double*, int, double*, double*, double*, double*, double*, double*);
    
	virtual double fx(lexer*, fdm_nhf*, double*, double*, int, double);
	virtual double fy(lexer*, fdm_nhf*, double*, double*, int, double);
	virtual double fz(lexer*, fdm_nhf*, double*, double*, int, double);
	void iqmin(lexer*, double*, double*, int);
	void jqmin(lexer*, double*, double*, int);
	void kqmin(lexer*, double*, double*, int);
	void iqmax(lexer*, double*, double*, int);
	void jqmax(lexer*, double*, double*, int);
	void kqmax(lexer*, double*, double*, int);

    
	double L,grad;
    double sig;

	double gradx, grady, gradz;
	double fu1,fv1,fw1,fu2,fv2,fw2;
    
    double ivel1,ivel2,jvel1,jvel2,kvel1,kvel2;
    double Pval,Qval;

    
    nhflow_flux *pflux;

};

#endif
