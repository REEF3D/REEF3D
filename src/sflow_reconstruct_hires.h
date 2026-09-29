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

#ifndef SFLOW_RECONSTRUCT_HIRES_H_
#define SFLOW_RECONSTRUCT_HIRES_H_

#include"sflow_reconstruct.h"
#include"increment.h"
#include"slice4.h"

class lexer;
class ghostcell;
class fdm2D;
class slice;
class patchBC_interface;

using namespace std;

// second-order MUSCL face reconstruction with slope limiter (A 211)
//   0: first-order (no slope), 1: van Leer, 2: Superbee, 3: van Albada

class sflow_reconstruct_hires final : public sflow_reconstruct, public increment
{
public:
	sflow_reconstruct_hires(lexer*,patchBC_interface*);
	virtual ~sflow_reconstruct_hires();

    void reconstruct_x(lexer*,ghostcell*,fdm2D*,slice&,slice&,slice&) override final;
    void reconstruct_y(lexer*,ghostcell*,fdm2D*,slice&,slice&,slice&) override final;
    void reconstruct_WL(lexer*,ghostcell*,fdm2D*) override final;

private:
    double limiter(lexer*, double, double);
    
    slice4 dfdx,dfdy;
    
    double dfdx_plus,dfdx_min,dfdy_plus,dfdy_min;
    double val,denom,r,phi;
    
    patchBC_interface *pBC;
};

#endif
