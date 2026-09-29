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

#include"iowave.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void iowave::nhflow_precalc_dirichlet(lexer *p, fdm_nhf *d, ghostcell *pgc)
{  
        double etaval=0.0;
        double Dact,uwave,vwave,wwave;
        
        p->wavetime = p->simtime;
        
        for(n=0;n<p->gcslin_count;n++)
        {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];
        
        xg=xgen(p);
        yg=ygen(p);
        
        etaval = wave_eta(p,pgc,xg,yg);
        
        eta(i-1,j) =  etaval;
        eta(i-2,j) =  etaval;
        eta(i-3,j) =  etaval;
        }
        
        count=0;
        for(n=0;n<p->gcin_count;n++)
		{
		i=p->gcin[n][0];
		j=p->gcin[n][1];
		k=p->gcin[n][2];

        
        x=xgen(p);
        y=ygen(p);
            
        if(p->A515==4)
        {
        // mass-conservative wavemaker:
        // - paddle kinematics sampled at the still-water sigma positions, so the
        //   layer fluxes integrate to the paddle's still-water stroke flux
        //   (piston: u*h, flap: theta_dot*int(z+l)dz over the still-water depth)
        // - wave flux over the still-water depth -> zero-mean flux, no net
        //   second-order (Stokes) mass input <u_p*eta> through the boundary
        // - current flux over the actual depth
        // - ghost velocity consistent with the actual depth: U = UH/D
        Dact = MAX(d->eta(i,j) + d->depth(i,j), p->A544);
        
        z = p->ZP[KP]*d->depth(i,j) + d->bed(i,j) - p->phimean;
        
        uwave = wave_u(p,pgc,x,y,z);
        vwave = wave_v(p,pgc,x,y,z);
        wwave = wave_w(p,pgc,x,y,z);
        
        // U
        UHval[count] = d->depth(i,j)*uwave + Dact*p->Ui;
        uval[count]  = UHval[count]/Dact;
        
        // V
        VHval[count] = d->depth(i,j)*vwave;
        vval[count]  = VHval[count]/Dact;
        
        // W
        wval[count]  = wwave;
        WHval[count] = Dact*wval[count];
        }
        
        else
        {
            if(p->A515==1)
            etaval = 0.0;
            
            if(p->A515==2)
            etaval = d->eta(i,j);
            
            if(p->A515==3)
            etaval = eta(i-1,j);

        z = p->ZSP[IJK]-p->phimean;

        // U
        uval[count] = wave_u(p,pgc,x,y,z) + p->Ui;
        UHval[count] = (etaval + d->depth(i,j))*uval[count];
        
        // V
        vval[count] = wave_v(p,pgc,x,y,z);
        VHval[count] = (etaval + d->depth(i,j))*vval[count];
        
        // W
        wval[count] = wave_w(p,pgc,x,y,z);
        WHval[count] = (etaval + d->depth(i,j))*wval[count];
        }

        ++count;
        }
}




