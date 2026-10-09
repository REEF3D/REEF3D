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

// grid quantities from the parcels: solid fraction, solid velocity, particle stress and its gradient
void CPM::grid_update(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s, double *PX, double *PY, double *PZ, double *PU, double *PV, double *PW)
{
    volfrac_update(p,pgc,s,PX,PY,PZ,PU,PV,PW);
    
    if(p->Q12==1)
    stress_snider(p,pgc,s);
    
    else if(p->Q12==2)
    stress_packedbed(p,pgc,s);
    
    else
    {
        BASELOOP
        Tau(i,j,k)=0.0;
        
        pgc->start4a(p,Tau,1);
        
        cmax=0.0;
    }
    
    stress_gradient(p,a,pgc,s);
    
    dispersion_update(p,a,pgc);
    
    exposure_update(p,a);
    
    // bed shear stress per column (bedload layer Q 58, log)
    if(p->S10!=2)
    bedload_columns(p,a,pgc,s);

    if(p->Q58==0)
    bagnold_update(p,a,pgc,s);
}

// adaptive sub-steps over the fluid time step:
// particle CFL number S 14 and wave speed of the particle stress Q 18, at most Q 28 sub-steps
double CPM::substep_size(lexer *p, ghostcell *pgc, double trem, int qs)
{
    double vmax=0.0;
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE)
    vmax = MAX(vmax, sqrt(P.U[n]*P.U[n] + P.V[n]*P.V[n] + P.W[n]*P.W[n]));
    
    vmax = pgc->globalmax(vmax);
    
    double dtlim = trem;
    
    if(vmax>1.0e-15)
    dtlim = MIN(dtlim, p->S14*hmin/vmax);
    
    if(cmax>1.0e-15)
    dtlim = MIN(dtlim, p->Q18*hmin/cmax);
    
    // random displacement: rms step sqrt(2 K dt) below half a cell
    if(p->Q52==1 && Ktmax>1.0e-15)
    dtlim = MIN(dtlim, 0.125*hmin*hmin/Ktmax);
    
    // last allowed sub-step
    if(qs>=p->Q28-1)
    {
        if(dtlim<trem && p->mpirank==0)
        cout<<"CPM warning: maximum number of sub-steps Q 28 reached, stable step "<<dtlim<<" < "<<trem<<endl;
        
        dtlim = trem;
    }
    
    // avoid a tiny last sub-step
    if(trem-dtlim < 0.1*dtlim)
    dtlim = trem;
    
    return dtlim;
}

void CPM::mppic_EE1(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s, turbulence *pturb)
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
        
        substep_euler(p,a,pgc,s,pturb,dts);
        
        trem -= dts;
        ++qs;
    }
    
    nsub = qs;
    dtsub = p->dt/double(MAX(qs,1));
}

// explicit Euler, point implicit drag and friction, grid quantities from grid_update
// the tentative step (XRK1, URK1) is checked against the walls and the free volume of the cells
void CPM::substep_euler(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s, turbulence *pturb, double dt)
{
    double fac;
    
    const bool layer = p->Q58>0 && p->S10!=2;
    
    if(layer)
    bedload_exchange(p,a,pgc,s,dt);
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE)
    {
        // sub-grid bedload layer
        if(layer && P.Hop[n]>0.0)
        {
            bedload_move(p,a,s,n,dt);
            continue;
        }
        
        advec_mppic(p, a, P, s, pturb,
                    P.X, P.Y, P.Z, P.U, P.V, P.W,
                    F, G, H, dt);

        // velocity update
        fac = 1.0/(1.0 + dt*Dpx);
        
        P.URK1[n] = (P.U[n] + dt*F)*fac;
        P.VRK1[n] = (P.V[n] + dt*G)*fac;
        P.WRK1[n] = (P.W[n] + dt*H)*fac;
        
        // friction
        if(p->Q12==2 && p->Q13==1)
        friction(p,a,P.X[n],P.Y[n],P.Z[n],P.URK1[n],P.VRK1[n],P.WRK1[n],dt,fac);
        
        // jammed bed
        P.URK1[n] *= 1.0-Hjam;
        P.VRK1[n] *= 1.0-Hjam;
        P.WRK1[n] *= 1.0-Hjam;
        
        // tentative position
        P.XRK1[n] = P.X[n] + dt*P.URK1[n];
        P.YRK1[n] = P.Y[n] + dt*P.VRK1[n];
        P.ZRK1[n] = P.Z[n] + dt*P.WRK1[n];
        
        // turbulent dispersion
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

    // walls, then grid-limited step
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

    // parallel transfer
    P.xchange(p,pgc,bedch,2);
}
