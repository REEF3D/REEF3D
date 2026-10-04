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

// gradient with linear extrapolation into the ghost cells at physical boundaries
void CPM::gradient(lexer *p, ghostcell *pgc, field &f, field &gx, field &gy, field &gz)
{
    double fm,fp;
    
    BASELOOP
    {
        // x
        fm = wallcell(p,i-1,j,k) ? 2.0*f(i,j,k)-f(i+1,j,k) : f(i-1,j,k);
        fp = wallcell(p,i+1,j,k) ? 2.0*f(i,j,k)-f(i-1,j,k) : f(i+1,j,k);
        gx(i,j,k) = (fp-fm)/(p->DXP[IM1]+p->DXP[IP]);
        
        // y
        if(p->j_dir==1)
        {
        fm = wallcell(p,i,j-1,k) ? 2.0*f(i,j,k)-f(i,j+1,k) : f(i,j-1,k);
        fp = wallcell(p,i,j+1,k) ? 2.0*f(i,j,k)-f(i,j-1,k) : f(i,j+1,k);
        gy(i,j,k) = (fp-fm)/(p->DYP[JM1]+p->DYP[JP]);
        }
        else
        gy(i,j,k) = 0.0;
        
        // z
        fm = wallcell(p,i,j,k-1) ? 2.0*f(i,j,k)-f(i,j,k+1) : f(i,j,k-1);
        fp = wallcell(p,i,j,k+1) ? 2.0*f(i,j,k)-f(i,j,k-1) : f(i,j,k+1);
        gz(i,j,k) = (fp-fm)/(p->DZP[KM1]+p->DZP[KP]);
    }
    
    pgc->start4a(p,gx,1);
    pgc->start4a(p,gy,1);
    pgc->start4a(p,gz,1);
}

// gradient of the particle normal stress
// at physical boundaries the stress is linearly extrapolated into the ghost cell,
// the wall carries the load of the bed
void CPM::stress_gradient(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s)
{
    gradient(p,pgc,Tau,dTx,dTy,dTz);
}

// fluid pressure gradient at the particles
// S 10 1: the bed is a solid boundary for the fluid, inside the bed the pore pressure is hydrostatic
void CPM::pressure_gradient(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s)
{
    double hs;
    
    BASELOOP
    {
        dPx(i,j,k) = (a->press(i+1,j,k) - a->press(i-1,j,k))/(p->DXP[IM1]+p->DXP[IP]);
        dPy(i,j,k) = p->j_dir==1 ? (a->press(i,j+1,k) - a->press(i,j-1,k))/(p->DYP[JM1]+p->DYP[JP]) : 0.0;
        dPz(i,j,k) = (a->press(i,j,k+1) - a->press(i,j,k-1))/(p->DZP[KM1]+p->DZP[KP]);
        
        if(p->S10!=2)
        {
            // inside the bed seen by the parcels: hydrostatic (with the bedload layer Q 58 the iso-surface)
            hs = heaviside(((p->Q58>0 && zsplit==0) ? Tiso(i,j,k) : a->topo(i,j,k)) - epsi);
            
            dPx(i,j,k) = hs*dPx(i,j,k) + (1.0-hs)*p->W1*p->W20;
            dPy(i,j,k) = hs*dPy(i,j,k) + (1.0-hs)*p->W1*p->W21;
            dPz(i,j,k) = hs*dPz(i,j,k) + (1.0-hs)*p->W1*p->W22;
        }
    }

    pgc->start4a(p,dPx,1);
    pgc->start4a(p,dPy,1);
    pgc->start4a(p,dPz,1);
}
