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
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"partres2.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include"turbulence.h"

void partres2::move_RK2(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s, turbulence *pturb)
{
    count_particles(p,a,pgc,s);
    
// RK step 1
    stress_tensor(p, pgc, s);
    
    ALOOP
    a->test(i,j,k) = Tau(i,j,k);
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE)
    {
        //cout<<"n: "<<n<<" P.index: "<<P.index<<endl;
        
        if(p->Q11==1)
        advec_plain(p, a, P, s, pturb, 
                        P.X, P.Y, P.Z, P.U, P.V, P.W,
                        F, G, H, 1.0);
        
        if(p->Q11==2)
        advec_pic(p, a, P, s, pturb, 
                        P.X, P.Y, P.Z, P.U, P.V, P.W,
                        F, G, H, 1.0);
                                         
        // Velocity update
        P.URK1[n] = P.U[n] + p->dtsed*F;
        P.VRK1[n] = P.V[n] + p->dtsed*G;
        P.WRK1[n] = P.W[n] + p->dtsed*H;
        
        // Position update
        P.XRK1[n] = P.X[n] + p->dtsed*P.URK1[n];
        P.YRK1[n] = P.Y[n] + p->dtsed*P.VRK1[n];
        P.ZRK1[n] = P.Z[n] + p->dtsed*P.WRK1[n];
    }
    // cellSum update
    cellSum_update(p,pgc,s,1);
    
    // parallel transfer
    P.xchange(p,pgc,1);
    
    
// RK step 2
    stress_tensor(p, pgc, s);
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE)
    {
        if(p->Q11==1)
        advec_plain(p, a, P, s, pturb, 
                        P.XRK1, P.YRK1, P.ZRK1, P.URK1, P.VRK1, P.WRK1,
                        F, G, H, 0.5);
                        
        if(p->Q11==2)
        advec_pic(p, a, P, s, pturb, 
                        P.XRK1, P.YRK1, P.ZRK1, P.URK1, P.VRK1, P.WRK1,
                        F, G, H, 0.5);
                        
        // Velocity update
        P.U[n] = 0.5*P.U[n] + 0.5*P.URK1[n] + 0.5*p->dtsed*F;
        P.V[n] = 0.5*P.V[n] + 0.5*P.VRK1[n] + 0.5*p->dtsed*G;
        P.W[n] = 0.5*P.W[n] + 0.5*P.URK1[n] + 0.5*p->dtsed*H;
        
        // Position update
        P.X[n] = 0.5*P.X[n] + 0.5*P.XRK1[n] + 0.5*p->dtsed*P.U[n];
        P.Y[n] = 0.5*P.Y[n] + 0.5*P.YRK1[n] + 0.5*p->dtsed*P.V[n];
        P.Z[n] = 0.5*P.Z[n] + 0.5*P.ZRK1[n] + 0.5*p->dtsed*P.W[n];
    }
    // cellSum update
    cellSum_update(p,pgc,s,2);
    
    // parallel transfer
    P.xchange(p, pgc, 2);
    
    //cellSum_full_update(p,pgc,s);
}