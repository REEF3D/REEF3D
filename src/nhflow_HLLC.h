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

#include"nhflow_convection.h"
#include"slice1.h"
#include"slice2.h"
#include"increment.h"

class nhflow_flux_build;

class patchBC_interface;
class ghostcell;

#ifndef NHFLOW_HLLC_H_
#define NHFLOW_HLLC_H_

using namespace std;

class nhflow_HLLC : public nhflow_convection, public increment
{

public:

	nhflow_HLLC (lexer*,ghostcell*,patchBC_interface*);
	virtual ~nhflow_HLLC();

    virtual void start(lexer*&, fdm_nhf*&, int, slice&);
    virtual void precalc(lexer*, fdm_nhf*, int, slice&);

private:

    void aij_U(lexer*, fdm_nhf*, int);
    void aij_V(lexer*, fdm_nhf*, int);
    void aij_W(lexer*, fdm_nhf*, int);
    void aij_E(lexer*, fdm_nhf*, int);
    
    void HLLC(lexer*, fdm_nhf*, double*, double*, double*, double*, double*, double*, double*, double*);
    void HLLC_E(lexer*, fdm_nhf*);
    
    void HLL(lexer*, fdm_nhf*&, double*, double*, double*, double*);
    
	double dx,dy,dz;
	double udir,vdir,wdir;
	double L;
    double denom;

    
    double FsS,FnS,FeS,FwS;

    ghostcell *pgc;
    patchBC_interface *pBC;
    nhflow_flux_build *pflux;
};

#endif
