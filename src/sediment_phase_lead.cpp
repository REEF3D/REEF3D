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


#include"sediment_f.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include"slice1.h"
#include"slice2.h"
#include"slice4.h"

/*--------------------------------------------------------------------
S 37 2 (NHFLOW): phase-lead near-bed velocity of Nielsen (2006) for intra-wave sediment transport,
    u_e = cos(phi) u + sin(phi)/omega du/dt,   phi = S 140 [deg], omega = 2 pi/T,
T = S 141, else S 58, else the period of the wave generation. Nielsen (1992, 2006) and Nielsen & Callaghan
(2003): the bed shear stress of sheet flow leads the free-stream velocity, so that waves with acceleration
skewness (steep fronts, saw-tooth waves in the surf zone) carry sediment onshore even with velocity
skewness zero (Hoefel & Elgar 2003, Duck94).
u is the velocity of the bottom sigma layer; the bed shear stress of S 16 1 (s->ub, s->vb) and the face
velocities that give the direction of the bed load (s->P, s->Q) use u_e. du/dt is taken between sediment
steps and low-pass filtered with the time constant T/50 (the non-hydrostatic pressure gives grid-scale
noise of the acceleration at the time-step scale) and limited to |du/dt| <= 0.3 g; no phase lead in water
shallower than max(0.1 m, 10 A 544) (swash: the velocities jump when cells wet or dry).
Called after fill_PQ_nhflow (prep_nhflow); Cartesian grids only (the curvilinear and AMR paths keep u).
--------------------------------------------------------------------*/

void sediment_f::phase_lead_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    double T = p->S141>0.0 ? p->S141 : (p->S58>0.0 ? p->S58 : (p->wTp>0.0 ? p->wTp : p->wT));
    
    k=0;
    
    // cell velocities of the bottom layer
    SLICELOOP4
    {
    s->ub(i,j) = d->U[IJK];
    s->vb(i,j) = d->V[IJK];
    }
    pgc->gcsl_start4(p,s->ub,1);
    pgc->gcsl_start4(p,s->vb,1);
    
    if(!(T>0.0))
    {
        if(pl_t0<-0.5 && p->mpirank==0)
        cout<<"S 37 2: no wave period (S 141, S 58 or the wave generation): no phase lead"<<endl;
        pl_t0 = 0.0;
        return;
    }
    
    if(pl_P0==nullptr)
    {
    pl_P0 = new slice1(p); pl_aP = new slice1(p);
    pl_Q0 = new slice2(p); pl_aQ = new slice2(p);
    pl_U0 = new slice4(p); pl_aU = new slice4(p);
    pl_V0 = new slice4(p); pl_aV = new slice4(p);
    
        for(i=p->imin; i<p->imin+p->imax; ++i)
        for(j=p->jmin; j<p->jmin+p->jmax; ++j)
        {
        (*pl_P0)(i,j) = s->P(i,j);  (*pl_aP)(i,j) = 0.0;
        (*pl_Q0)(i,j) = s->Q(i,j);  (*pl_aQ)(i,j) = 0.0;
        (*pl_U0)(i,j) = s->ub(i,j); (*pl_aU)(i,j) = 0.0;
        (*pl_V0)(i,j) = s->vb(i,j); (*pl_aV)(i,j) = 0.0;
        }
    
    pl_t0 = p->simtime;
    
        if(p->mpirank==0)
        cout<<"S 37 2: phase-lead near-bed velocity (Nielsen 2006), phi = "<<p->S140<<" deg, T = "<<T<<" s"<<endl;
    }
    
    const double dts = p->simtime - pl_t0;
    pl_t0 = p->simtime;
    
    const double phi = p->S140*(PI/180.0);
    const double cp = cos(phi);
    const double sp = sin(phi)*T/(2.0*PI);
    const double w = dts>0.0 ? MIN(1.0, dts/(T/50.0)) : 0.0;
    const double idt = dts>1.0e-12 ? 1.0/dts : 0.0;
    
    // limits: no phase lead in thin water (swash films, wetting and drying: the face and cell velocities jump
    // when a cell wets or dries, and du/dt with them) and |du/dt| <= 0.3 g
    const double hpl = MAX(0.1, 10.0*p->A544);
    const double amax = 0.3*fabs(p->W22);
    const int ie = p->imin+p->imax, je = p->jmin+p->jmax;
    
    auto lead = [&](double &u, double &u0, double &af, bool wet)
    {
        double a = wet ? (u - u0)*idt : 0.0;
        a = MAX(-amax, MIN(amax, a));
        af = wet ? af + w*(a - af) : 0.0;
        u0 = u;
        if(wet)
        u = cp*u + sp*af;
    };
    
    // the indices of all slice points incl. the ghost points (as transport_direction): no exchange needed
    for(i=p->imin; i<ie; ++i)
    for(j=p->jmin; j<je; ++j)
    {
    const double hc = s->waterlevel(i,j);
    const double hx = i+1<ie ? MIN(hc, s->waterlevel(i+1,j)) : hc;
    const double hy = j+1<je ? MIN(hc, s->waterlevel(i,j+1)) : hc;
    
    lead(s->P(i,j),  (*pl_P0)(i,j), (*pl_aP)(i,j), hx>hpl);
    lead(s->Q(i,j),  (*pl_Q0)(i,j), (*pl_aQ)(i,j), hy>hpl);
    lead(s->ub(i,j), (*pl_U0)(i,j), (*pl_aU)(i,j), hc>hpl);
    lead(s->vb(i,j), (*pl_V0)(i,j), (*pl_aV)(i,j), hc>hpl);
    }
}
