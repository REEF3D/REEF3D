/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#ifndef SFLOW_HLL_H_
#define SFLOW_HLL_H_

#include"sflow_convection.h"
#include"slice1.h"
#include"slice2.h"
#include"increment.h"

class sflow_flux_build;

class patchBC_interface;
class ghostcell;

using namespace std;

class sflow_HLL : public sflow_convection, public increment
{

public:

	sflow_HLL (lexer*,ghostcell*,patchBC_interface*);
	virtual ~sflow_HLL();

    virtual void start(lexer*&, fdm2D*&, int, slice&);
    virtual void precalc(lexer*, fdm2D*, int, slice&);

private:

    void aij_U(lexer*&, fdm2D*&, int);
    void aij_V(lexer*&, fdm2D*&, int);
    void aij_W(lexer*&, fdm2D*&, int);
    void aij_E(lexer*&, fdm2D*&, int);
    
    void HLL(lexer*&, fdm2D*&, slice &, slice &, slice &, slice &);
    void HLL_E(lexer*&, fdm2D*&);
    
	double dx,dy,dz;
	double udir,vdir,wdir;
	double L;
    double denom;

    double ivel1,ivel2,jvel1,jvel2,kvel1,kvel2;

    ghostcell *pgc;
    patchBC_interface *pBC;
    sflow_flux_build *pflux;
};

#endif
