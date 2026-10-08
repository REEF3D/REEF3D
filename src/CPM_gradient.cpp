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
// S 10 1: the bed is a solid boundary for the fluid, inside the bed the pore pressure is hydrostatic,
// or with the seepage flow Q 69 the pore pressure of the Darcy flow (CPM_seepage.cpp).
// The fluid solves no pressure in the cells of its bed (topo < 0) and of solid bodies: the values there
// are what the initialisation and the ghost-cell updates left, not a pore pressure. Next to such a cell
// the gradient takes the pore pressure of the cell instead (hydrostatic continuation of the cell itself,
// or the seepage pressure), and inside it the gradient of the pore pressure; otherwise the parcels near a
// sloping bed surface feel a lateral pressure gradient that depends on the initial pressure field (with
// the initial pressure of the fluid the bed was held by a suction, a 30 degree wedge stood only with it).
void CPM::pressure_gradient(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s)
{
    double hs;
    const bool seep = p->Q69>0 && p->S10!=2;
    
    if(seep)
    seepage_update(p,a,pgc);
    
    auto nofluid = [&](int ii, int jj, int kk) -> bool
    {
        return a->topo(ii,jj,kk)<0.0 || (p->solidread>0 && a->solid(ii,jj,kk)<0.0);
    };
    
    // pressure of the neighbour (ii,jj,kk) at the distance (dx,dy,dz) from the cell (i,j,k)
    auto pnb = [&](int ii, int jj, int kk, double dx, double dy, double dz) -> double
    {
        if(p->S10!=2 && nofluid(ii,jj,kk))
        {
            if(seep && seep_bed(p,a,ii,jj,kk))
            return Pse(ii,jj,kk) + p->W1*(p->W20*p->XP[ii+marge] + p->W21*p->YP[jj+marge] + p->W22*p->ZP[kk+marge]);
            
            return a->press(i,j,k) + p->W1*(p->W20*dx + p->W21*dy + p->W22*dz);
        }
        
        return a->press(ii,jj,kk);
    };
    
    // seepage: gradient of p* in a cell of the bed; walls and solid bodies mirror the cell
    auto pse = [&](int ii, int jj, int kk) -> double
    {
        if(wallcell(p,ii,jj,kk) || (p->solidread>0 && a->solid(ii,jj,kk)<0.0))
        return Pse(i,j,k);
        
        return Pse(ii,jj,kk);
    };
    
    BASELOOP
    {
        if(p->S10!=2 && nofluid(i,j,k))
        {
            dPx(i,j,k) = p->W1*p->W20;
            dPy(i,j,k) = p->W1*p->W21;
            dPz(i,j,k) = p->W1*p->W22;
            
            // excess pore pressure of the seepage: the force on the bed per unit volume, taken by the grains
            if(seep && seep_bed(p,a,i,j,k))
            {
                double tdiv = MAX(Ts(i,j,k), 0.5*theta_bed);
                
                dPx(i,j,k) += (pse(i+1,j,k) - pse(i-1,j,k))/(p->DXP[IM1]+p->DXP[IP])/tdiv;
                
                if(p->j_dir==1)
                dPy(i,j,k) += (pse(i,j+1,k) - pse(i,j-1,k))/(p->DYP[JM1]+p->DYP[JP])/tdiv;
                
                dPz(i,j,k) += Gsz(i,j,k)/tdiv;
            }
            
            continue;
        }
        
        dPx(i,j,k) = (pnb(i+1,j,k,p->DXP[IP],0.0,0.0) - pnb(i-1,j,k,-p->DXP[IM1],0.0,0.0))/(p->DXP[IM1]+p->DXP[IP]);
        dPy(i,j,k) = p->j_dir==1 ? (pnb(i,j+1,k,0.0,p->DYP[JP],0.0) - pnb(i,j-1,k,0.0,-p->DYP[JM1],0.0))/(p->DYP[JM1]+p->DYP[JP]) : 0.0;
        dPz(i,j,k) = (pnb(i,j,k+1,0.0,0.0,p->DZP[KP]) - pnb(i,j,k-1,0.0,0.0,-p->DZP[KM1]))/(p->DZP[KM1]+p->DZP[KP]);
        
        if(p->S10!=2)
        {
            // inside the bed seen by the parcels (with the bedload layer Q 58 the iso-surface): the pore pressure,
            // switched at the bed surface (a smooth step over the interface width of the free surface,
            // F 45 h, would replace the dynamic pressure gradient of the flow in the first cells above the bed)
            hs = ((p->Q58>0 && zsplit==0) ? Tiso(i,j,k) : a->topo(i,j,k)) >= 0.0 ? 1.0 : 0.0;
            
            dPx(i,j,k) = hs*dPx(i,j,k) + (1.0-hs)*p->W1*p->W20;
            dPy(i,j,k) = hs*dPy(i,j,k) + (1.0-hs)*p->W1*p->W21;
            dPz(i,j,k) = hs*dPz(i,j,k) + (1.0-hs)*p->W1*p->W22;
        }
    }

    pgc->start4a(p,dPx,1);
    pgc->start4a(p,dPy,1);
    pgc->start4a(p,dPz,1);
}
