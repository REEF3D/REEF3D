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

void CPM::mppic_RK2(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s, turbulence *pturb)
{
    count_particles(p,a,pgc,s);

    pressure_gradient(p,a,pgc,s);
    
    double trem = p->dt;
    double dts;
    int qs=0;
    
    while(trem>1.0e-10*p->dt)
    {
        grid_update(p,a,pgc,s,P.X,P.Y,P.Z,P.U,P.V,P.W);
        
        dts = substep_size(p,pgc,trem,qs);
        
        substep_rk2(p,a,pgc,s,pturb,dts);
        
        trem -= dts;
        ++qs;
    }
    
    nsub = qs;
    dtsub = p->dt/double(MAX(qs,1));
}

// Heun / SSP-RK2, point implicit drag and friction in each stage
void CPM::substep_rk2(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s, turbulence *pturb, double dt)
{
    double fac;
    
    const bool layer = p->Q58>0 && p->S10!=2;
    
    if(layer)
    bedload_exchange(p,a,pgc,s,dt);
    
 // ------------------------
    // RK step 1, grid quantities from grid_update
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE)
    {
        // sub-grid bedload layer: one step for both stages
        if(layer && P.Hop[n]>0.0)
        {
            bedload_move(p,a,s,n,dt);
            continue;
        }
        
        advec_mppic(p, a, P, s, pturb,
                    P.X, P.Y, P.Z, P.U, P.V, P.W,
                    F, G, H, dt);

        fac = 1.0/(1.0 + dt*Dpx);
        
        P.URK1[n] = (P.U[n] + dt*F)*fac;
        P.VRK1[n] = (P.V[n] + dt*G)*fac;
        P.WRK1[n] = (P.W[n] + dt*H)*fac;
        
        if(p->Q12==2 && p->Q13==1)
        friction(p,a,P.X[n],P.Y[n],P.Z[n],P.URK1[n],P.VRK1[n],P.WRK1[n],dt,fac);
        
        // jammed bed
        P.URK1[n] *= 1.0-Hjam;
        P.VRK1[n] *= 1.0-Hjam;
        P.WRK1[n] *= 1.0-Hjam;
        
        P.XRK1[n] = P.X[n] + dt*P.URK1[n];
        P.YRK1[n] = P.Y[n] + dt*P.VRK1[n];
        P.ZRK1[n] = P.Z[n] + dt*P.WRK1[n];
    }
    
    boundcheck(p,1);
    
    if(p->Q19==1)
    limiter(p,a,pgc,P.X,P.Y,P.Z,P.XRK1,P.YRK1,P.ZRK1,P.URK1,P.VRK1,P.WRK1);
    
    periodic_wrap(p,P.X,P.Y,P.XRK1,P.YRK1);

    P.xchange(p,pgc,bedch,1);

    // ------------------------
    // RK step 2
    grid_update(p,a,pgc,s,P.XRK1,P.YRK1,P.ZRK1,P.URK1,P.VRK1,P.WRK1);

    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE)
    {
        if(layer && P.Test[n]<0.0)
        continue;
        
        advec_mppic(p, a, P, s, pturb,
                    P.XRK1, P.YRK1, P.ZRK1, P.URK1, P.VRK1, P.WRK1,
                    F, G, H, dt);

        fac = 1.0/(1.0 + 0.5*dt*Dpx);
        
        P.U[n] = (0.5*P.U[n] + 0.5*P.URK1[n] + 0.5*dt*F)*fac;
        P.V[n] = (0.5*P.V[n] + 0.5*P.VRK1[n] + 0.5*dt*G)*fac;
        P.W[n] = (0.5*P.W[n] + 0.5*P.WRK1[n] + 0.5*dt*H)*fac;
        
        if(p->Q12==2 && p->Q13==1)
        friction(p,a,P.XRK1[n],P.YRK1[n],P.ZRK1[n],P.U[n],P.V[n],P.W[n],0.5*dt,fac);
        
        // jammed bed
        P.U[n] *= 1.0-Hjam;
        P.V[n] *= 1.0-Hjam;
        P.W[n] *= 1.0-Hjam;
        
        // tentative position in XRK1, velocity in URK1
        P.XRK1[n] = 0.5*P.X[n] + 0.5*P.XRK1[n] + 0.5*dt*P.U[n];
        P.YRK1[n] = 0.5*P.Y[n] + 0.5*P.YRK1[n] + 0.5*dt*P.V[n];
        P.ZRK1[n] = 0.5*P.Z[n] + 0.5*P.ZRK1[n] + 0.5*dt*P.W[n];
        
        P.URK1[n] = P.U[n];
        P.VRK1[n] = P.V[n];
        P.WRK1[n] = P.W[n];
        
        // turbulent dispersion, once per step from the start position
        if(p->Q52==1 && !(layer && bedload_nodisp(p,a,n)))
        {
            double ddx,ddy,ddz;
            dispersion(p,P.X[n],P.Y[n],P.Z[n],ddx,ddy,ddz,dt);
            P.XRK1[n] += ddx;
            P.YRK1[n] += ddy;
            P.ZRK1[n] += ddz;
        }
        
        // sub-grid bedload layer: a parcel settling from the flow onto the bed is deposited
        if(layer)
        bedload_settle(p,a,s,n);
    }

    boundcheck(p,1);
    
    if(p->Q19==1)
    limiter(p,a,pgc,P.X,P.Y,P.Z,P.XRK1,P.YRK1,P.ZRK1,P.URK1,P.VRK1,P.WRK1);
    
    periodic_wrap(p,P.X,P.Y,P.XRK1,P.YRK1);
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE)
    {
        P.X[n] = P.XRK1[n];
        P.Y[n] = P.YRK1[n];
        P.Z[n] = P.ZRK1[n];
        P.U[n] = P.URK1[n];
        P.V[n] = P.VRK1[n];
        P.W[n] = P.WRK1[n];
    }

    P.xchange(p,pgc,bedch,2);
}
