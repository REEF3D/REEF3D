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


#include"CPM.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

/*--------------------------------------------------------------------
turbulent dispersion of the parcels (Q 52 1): random displacement model

    dx = grad(K) dt + sqrt(2 K dt) xi,     xi ~ N(0,1) per direction

K = nu_t/Sc_t (Q 53) is the eddy diffusivity of the fluid. The drift grad(K) keeps
a well-mixed suspension well mixed in inhomogeneous turbulence, so the ensemble of the
parcels follows d c/dt + div(c u_p) = div(K grad c): with the settling velocity this gives
the Rouse profile. The parcel velocity is not changed, the displacement is added to the
tentative position before the grid-limited step, so the packing limit holds.

K is taken in the water only (phi > 0) and is damped towards the bed, it vanishes where
the contact network begins to carry the grains (the gap 0.5 theta_bed of the overburden
stress), so the grains of the bed and the bedload layer are not kicked by the random walk:

    K = nu_t/Sc_t  max(0, 1 - theta/(0.5 theta_bed))
--------------------------------------------------------------------*/

void CPM::dispersion_update(lexer *p, fdm *a, ghostcell *pgc)
{
    if(p->Q52!=1)
    return;
    
    double Sc = p->Q53>0.0 ? p->Q53 : 1.0;
    double kmax=0.0;
    
    BASELOOP
    {
        double damp = MAX(0.0, MIN(1.0, 1.0 - Ts(i,j,k)/MAX(0.5*theta_bed,1.0e-6)));
        double water = a->phi(i,j,k)>0.0 ? 1.0 : 0.0;
        
        Kt(i,j,k) = water*damp*MAX(a->eddyv(i,j,k),0.0)/Sc;
        
        if(p->flag4[IJK]<0)
        Kt(i,j,k) = 0.0;
        
        kmax = MAX(kmax,Kt(i,j,k));
    }
    
    pgc->start4a(p,Kt,1);
    
    gradient(p,pgc,Kt,dKx,dKy,dKz);
    
    Ktmax = pgc->globalmax(kmax);
}

void CPM::dispersion(lexer *p, double xp, double yp, double zp, double &dx, double &dy, double &dz, double dt)
{
    dx=dy=dz=0.0;
    
    if(p->Q52!=1)
    return;
    
    double K = MAX(0.0, p->ccipol4a(Kt,xp,yp,zp));
    
    if(K<=0.0)
    return;
    
    double sd = sqrt(2.0*K*dt);
    
    dx = p->ccipol4a(dKx,xp,yp,zp)*dt + sd*gauss(rng);
    dz = p->ccipol4a(dKz,xp,yp,zp)*dt + sd*gauss(rng);
    
    if(p->j_dir==1)
    dy = p->ccipol4a(dKy,xp,yp,zp)*dt + sd*gauss(rng);
}
