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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#include "partres.h"
#include "particles_obj.h"
#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"

void partres::move_RK2(lexer *p, fdm &a, ghostcell &pgc, particles_obj &PP, sediment_fdm &s, turbulence &pturb)
{
    double du,dv,dw;
    
    // RK step1 
    for(size_t n=0;n<PP.loopindex;n++)
    if(PP.Flag[n]>0)
    {
        PP.URK1[n] = PP.U[n] + p->dtsed*du;
        PP.VRK1[n] = PP.V[n] + p->dtsed*dv;
        PP.WRK1[n] = PP.W[n] + p->dtsed*dw;
        
        // Pos update
        PP.X[n] += p->dtsed*PP.U[n];
        PP.Y[n] += p->dtsed*PP.V[n];
        PP.Z[n] += p->dtsed*PP.W[n];

        // Particel sum update
        cellSum[IJK]-=PP.PackingFactor[n];
        i=p->posc_i(PP.X[n]);
        j=p->posc_j(PP.Y[n]);
        k=p->posc_k(PP.Z[n]);
        cellSum[IJK]+=PP.PackingFactor[n];
    }
    
}