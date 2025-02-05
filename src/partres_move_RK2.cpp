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
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"partres.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include"turbulence.h"

void partres::move_RK2(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s, turbulence *pturb)
{
    count_particles(p,a,pgc,s);
    pressure_gradient(p,a,pgc,s);
    
// RK step 1
    stress_tensor(p,pgc,s);
    stress_gradient(p,a,pgc,s);
    
    ALOOP
    a->test(i,j,k) = cellSum(i,j,k)/P.ParcelFactor;
    //a->test(i,j,k) = Ts(i,j,k);
   // a->test(i,j,k) = rf(p,p->pos_x(),p->pos_y());
    //a->test(i,j,k) = (Tau(i,j,k+1) - Tau(i,j,k-1))/(p->DZP[KM1]+p->DZP[KP]);
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE)
    {        
        if(p->Q11==1)
        advec_plain(p, a, P, s, pturb, 
                        P.X, P.Y, P.Z, P.U, P.V, P.W,
                        F, G, H, 1.0);
        
        if(p->Q11==2)
        advec_mppic(p, a, P, s, pturb, 
                        P.X, P.Y, P.Z, P.U, P.V, P.W,
                        F, G, H, 1.0);
                                         
        // Velocity update
        P.URK1[n] = P.U[n] + p->dtsed*F;
        P.VRK1[n] = P.V[n] + p->dtsed*G;
        P.WRK1[n] = P.W[n] + p->dtsed*H;
        
        // Position update
        P.XRK1[n] = P.X[n] + p->dtsed*P.URK1[n];
        P.YRK1[n] = P.Y[n] + p->dtsed*P.VRK1[n];
        //P.ZRK1[n] = P.Z[n] + p->dtsed*P.WRK1[n];
        P.ZRK1[n] = MAX(MIN(P.Z[n] + p->dtsed*P.WRK1[n],0.49),-0.3);
    }
    
    
    // cellSum update
    cellSum_full_update(p,pgc,s,1);
    
    boundcheck(p,a,pgc,s,1);
    bedchange_update(p,pgc,s,1);
    bedchange(p,a,pgc,s,1);
    
    // parallel transfer
    P.xchange(p,pgc,bedch,1);

// RK step 2
    stress_tensor(p, pgc, s);
    stress_gradient(p,a,pgc,s);
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE)
    {
        if(p->Q11==1)
        advec_plain(p, a, P, s, pturb, 
                        P.XRK1, P.YRK1, P.ZRK1, P.URK1, P.VRK1, P.WRK1,
                        F, G, H, 0.5);
                        
        if(p->Q11==2)
        advec_mppic(p, a, P, s, pturb, 
                        P.XRK1, P.YRK1, P.ZRK1, P.URK1, P.VRK1, P.WRK1,
                        F, G, H, 0.5);
                        
        // Velocity update
        P.U[n] = 0.5*P.U[n] + 0.5*P.URK1[n] + 0.5*p->dtsed*F;
        P.V[n] = 0.5*P.V[n] + 0.5*P.VRK1[n] + 0.5*p->dtsed*G;
        P.W[n] = 0.5*P.W[n] + 0.5*P.WRK1[n] + 0.5*p->dtsed*H;
        
        // Position update
        P.X[n] = 0.5*P.X[n] + 0.5*P.XRK1[n] + 0.5*p->dtsed*P.U[n];
        P.Y[n] = 0.5*P.Y[n] + 0.5*P.YRK1[n] + 0.5*p->dtsed*P.V[n];
        //P.Z[n] = 0.5*P.Z[n] + 0.5*P.ZRK1[n] + 0.5*p->dtsed*P.W[n];
        P.Z[n] = MAX(MIN(0.5*P.Z[n] + 0.5*P.ZRK1[n] + 0.5*p->dtsed*P.W[n],0.49),-0.3);
    }
    
    
    // cellSum update
    cellSum_full_update(p,pgc,s,2);
    
    boundcheck(p,a,pgc,s,2);
    bedchange_update(p,pgc,s,2);
    bedchange(p,a,pgc,s,2);
    
    // parallel transfer
    P.xchange(p, pgc,bedch,2);
    
    
    
}