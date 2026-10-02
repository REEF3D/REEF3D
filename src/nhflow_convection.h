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

#ifndef NHFLOW_CONVECTION_H_
#define NHFLOW_CONVECTION_H_

class lexer;
class fdm_nhf;
class slice;

using namespace std;

// mesh refinement (nhflow_amr): called with the face fluxes of variable ipol (1 UH, 2 VH, 3 WH:
// Fx, Fy; 4 continuity: FEx, FEy) after their ghost cells are set and before the divergence
class nhflow_flux_hook
{
public:
    virtual void flux_hook(lexer*, fdm_nhf*, int, int, double*, double*)=0;
};

class nhflow_convection
{
public:

    virtual void start(lexer*&, fdm_nhf*&, int, slice&, double*)=0;
    virtual void precalc(lexer*, fdm_nhf*, int, slice&)=0;
    
    void set_hook(nhflow_flux_hook *h, int id) { phook=h; hook_id=id; }
    nhflow_flux_hook *phook = nullptr;
    int hook_id = -1;

};

#endif



