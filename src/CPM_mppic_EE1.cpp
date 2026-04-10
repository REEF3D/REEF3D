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
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"CPM.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include"turbulence.h"

void CPM::mppic_EE1(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s, turbulence *pturb)
{
    count_particles(p,a,pgc,s);
    
    press_lithostatic(p,a,pgc,s);
    
    pressure_gradient(p,a,pgc,s);
    
    ALOOP
    test(i,j,k) = dTz(i,j,k);
    
    // stress and cellSum update
    volfrac_update(p,pgc,s,P.X,P.Y,P.Z);
    stress_snider(p,pgc,s);
    stress_gradient(p,a,pgc,s);
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE)
    {

        advec_mppic_step1(p, a, P, s, pturb,
                    P.X, P.Y, P.Z, P.U, P.V, P.W,
                    F, G, H, 0.5);

        // Velocity update 1
        P.U[n] = (P.U[n] + p->dtsed*F)/(1.0 + p->dtsed*Dpx);
        P.V[n] = (P.V[n] + p->dtsed*G)/(1.0 + p->dtsed*Dpy);
        P.W[n] = (P.W[n] + p->dtsed*H)/(1.0 + p->dtsed*Dpz);
        
        // Position update
        P.X[n] = P.X[n] + p->dtsed*P.U[n];
        P.Y[n] = P.Y[n] + p->dtsed*P.V[n];
        P.Z[n] = P.Z[n] + p->dtsed*P.W[n];
    }
    /*
    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE)
    {
        // advec 2
        advec_mppic_step2(p, a, P, s, pturb,
                    P.X, P.Y, P.Z, P.U, P.V, P.W,
                    F, G, H, 0.5);
                    
        //F=G=H=0.0;

        // Velocity update 2
        P.U[n] += 0.5*p->dtsed*F;
        P.V[n] += 0.5*p->dtsed*G;
        P.W[n] += 0.5*p->dtsed*H;

        // Position update
        P.X[n] = 0.5*P.X[n] + 0.5*P.XRK1[n] + 0.5*p->dtsed*P.U[n];
        P.Y[n] = 0.5*P.Y[n] + 0.5*P.YRK1[n] + 0.5*p->dtsed*P.V[n];
        P.Z[n] = 0.5*P.Z[n] + 0.5*P.ZRK1[n] + 0.5*p->dtsed*P.W[n];
    }*/

    boundcheck(p,2);
    //bedchange_update(p,pgc,2);
    //bedchange(p,a,pgc,s,2);

    // parallel transfer
    P.xchange(p, pgc,bedch,2);
}
