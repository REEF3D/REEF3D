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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#ifndef NHFLOW_SOLID_FLUX_H_
#define NHFLOW_SOLID_FLUX_H_

// Impermeable solid cells for the continuity fluxes (d->solid_flux, set by the
// FEM coupling, Z 30, and by nhflow_forcing for direct forcing solids, A 521 0):
// the fluxes FEx, FEy of the free surface equation are zero at every face of a
// solid cell (p->DF < 0, or p->DF == -1, see below), layer by layer. The face
// states of the momentum fluxes are already wall states there
// (nhflow_fsf_f::wetdry_fluxes), but the dissipative term of the HLL flux, Sn Ss (Dn - Ds), still carries water
// into and through a solid column. Applied after the fluxes are built, before
// their halo exchange.

#include"lexer.h"
#include"fdm_nhf.h"
#include"increment.h"

inline void nhflow_solid_flux(lexer *p, fdm_nhf *d)
{
    int i,j,k;
    
    // 1: every cell with DF < 0 (FEM structures, Z 30)
    // 2: fixed direct forcing solids only (DF == -1: A 590 STL, A 581-589 primitives, A 580 grid solids),
    //    floating bodies (DF == -2) keep their continuity fluxes
    const int mode = d->solid_flux;
    auto sol = [p,mode](int q){return mode==2 ? p->DF[q]==-1 : p->DF[q]<0;};

    LOOP
    {
        if(sol(IJK) || sol(Ip1JK))
        d->FEx[IJK] = 0.0;

        if(sol(IJK) || sol(Im1JK))
        d->FEx[Im1JK] = 0.0;

        if(p->j_dir==1)
        {
            if(sol(IJK) || sol(IJp1K))
            d->FEy[IJK] = 0.0;

            if(sol(IJK) || sol(IJm1K))
            d->FEy[IJm1K] = 0.0;
        }
    }
}

#endif
