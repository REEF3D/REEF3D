/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
        
        p->wavetime = p->simtime + p->dt;
        
        for(n=0;n<p->gcslin_count;n++)
        {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];
        
        eta0(i,j)   =  eta1(i,j);
        eta0(i-1,j) =  eta1(i,j);
        eta0(i-2,j) =  eta1(i,j);
        eta0(i-3,j) =  eta1(i,j);
        
        xg=xgen(p);
        yg=ygen(p);
        
        eta1(i,j) = wave_eta(p,pgc,xg,yg);
        
        /*eta1(i-1,j) =  d->eta(i,j);
        eta1(i-2,j) =  d->eta(i,j);
        eta1(i-3,j) =  d->eta(i,j);*/
        
        eta1(i-1,j) =  eta1(i,j);
        eta1(i-2,j) =  eta1(i,j);
        eta1(i-3,j) =  eta1(i,j);
        }
        
        count=0;
        for(n=0;n<p->gcin_count;n++)
		{
		i=p->gcin[n][0];
		j=p->gcin[n][1];
		k=p->gcin[n][2];
        
        uval0[count] = uval1[count];
        vval0[count] = vval1[count];
        wval0[count] = wval1[count];
        
        UHval0[count] = UHval1[count];
        VHval0[count] = VHval1[count];
        WHval0[count] = WHval1[count];
        
        x=xgen(p);
        y=ygen(p);
            
        etaval = eta0(i,j);
        
        if(p->B92>=20 && p->B92<=29)
        etaval = 0.0;

        z = p->ZSP[IJK]-p->phimean;

        // U
        uval1[count] = wave_u(p,pgc,x,y,z) + p->Ui;
        UHval1[count] = (etaval + d->depth(i,j))*uval1[count];
        
        // V
        vval1[count] = wave_v(p,pgc,x,y,z);
        VHval1[count] = (etaval + d->depth(i,j))*vval1[count];
        
        // W
        wval1[count] = wave_w(p,pgc,x,y,z);
        VHval1[count] = (etaval + d->depth(i,j))*wval1[count];

        ++count;
        }
}




