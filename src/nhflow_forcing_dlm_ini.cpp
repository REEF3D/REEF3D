/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"nhflow_forcing.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void nhflow_forcing::dlm_forcing_ini(lexer *p, ghostcell *pgc)
{
    Ne = 10;
    Np = Ne + 1;

    
    p->Darray(EL_L, p->A584);
    p->Darray(EL_dx, p->A584);
    p->Darray(EL_X, p->A584, Np);
    p->Darray(EL_Y, p->A584, Np);
    p->Darray(EL_Z, p->A584, Np);
    p->Darray(EL_V, p->A584, Np);
    p->Iarray(EL_f, p->A584, Np);
    p->Darray(EL_FX, p->A584, Np);
    p->Darray(EL_FY, p->A584, Np);
    p->Darray(EL_FZ, p->A584, Np);
    
    // Initialise parameter
    for(int n=0; n<entity_count; ++n)
    for(int q=0; q<Np; ++q)
    EL_f[n][q] = 0;
    
    for(int n=0; n<p->A584; ++n)
    if(p->A584_xc[n]>=p->originx && p->A584_xc[n]<p->endx
    && p->A584_yc[n]>=p->originy && p->A584_yc[n]<p->endy)
    {
        EL_L[n] = p->A584_ze[n] - p->A584_zs[n];
        EL_dx[n] = (p->A584_ze[n] - p->A584_zs[n]) / double(Ne);
        
        for(int q=0; q<Np; ++q)
        {
        // if on local proc
        EL_X[n][q] = p->A584_xc[n];
        EL_Y[n][q] = p->A584_yc[n];
        EL_Z[n][q] = p->A584_zs[n] + double(n)*(p->A584_ze[n] - p->A584_zs[n])/double(Ne);
        EL_V[n][q] = PI*p->A584_r[n]*p->A584_r[n]*EL_dx[n];
        
        EL_f[n][q] = 1;
        }
    }
        
}