/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#ifndef WENO_FLUX_NUG_H_
#define WENO_FLUX_NUG_H_

#include"convection.h"
#include"weno_nug_func.h"

class flux;

using namespace std;

class weno_flux_nug final : public convection, public weno_nug_func
{
public:
	weno_flux_nug(lexer*);
	virtual ~weno_flux_nug();

	void start(lexer*,fdm*,field&,int,field&,field&,field&) override final;

    // the face value of the cell scalar b at the face i+1/2 (d 0), j+1/2 (1), k+1/2 (2) of cell
    // (i,j,k) for the face velocity vel (as start with ipol 4: the conservative flux vel x value;
    // the mesh refinement takes the fluxes through the faces at its patch edges from it)
    double face_value(lexer*, fdm*, field &b, field &vel, int d, int i, int j, int k, double vel_face);

private:
    double aij(lexer*, fdm*, field&, int,field&,field&,field&,double*,double*,double*);
    
	double fx(lexer*, fdm*, field&, field&, int, double);
	double fy(lexer*, fdm*, field&, field&, int, double);
	double fz(lexer*, fdm*, field&, field&, int, double);
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
