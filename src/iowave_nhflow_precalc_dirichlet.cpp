/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

void iowave::nhflow_precalc_dirichlet(lexer *p,fdm_nhf *d, ghostcell *pgc)
{       
        count=0;
		for(n=0;n<p->gcslin_count;n++)
        {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];
        
        xg=xgen(p);
        yg=ygen(p);
        x1=xgen1(p);
        y2=ygen2(p);
        

        eta(i,j) = wave_eta(p,pgc,xg,yg);
        etaval[count] = eta(i,j);
        ++count;
        }
        
        
        count=0;
		for(n=0;n<p->gcin_count;n++)
		{
		i=p->gcin[n][0];
		j=p->gcin[n][1];
		k=p->gcin[n][2];
        
        //cout<<"precalc: "<<i<<" "<<j<<" "<<k<<endl;
        x=xgen(p);
        y=ygen(p);
        x1=xgen1(p);
        y2=ygen2(p);
                
        z = p->ZSP[IJK]-p->phimean;

        // U
        if(z<=eta(i,j)+epsi)
        {
        uval[count] = wave_u(p,pgc,x1,y,z);
        UHval[count] = (eta(i,j) + d->depth(i,j))*(wave_u(p,pgc,xg,yg,z) + p->Ui);
        }
        
        if(z>eta(i,j)+epsi)
        uval[count] = 0.0;
        
        // V
        if(z<=eta(i,j)+epsi)
        {
        vval[count] = wave_v(p,pgc,x,y2,z);
        VHval[count] = (eta(i,j) + d->depth(i,j))*wave_v(p,pgc,xg,yg,z);
        }
            
        if(z>eta(i,j)+epsi)
        vval[count] = 0.0;
        
        // W
        if(z<=eta(i,j)+epsi)
        {
        wval[count] = wave_w(p,pgc,x,y,z);
        VHval[count] = (eta(i,j) + d->depth(i,j))*wave_v(p,pgc,xg,yg,z);
        }
            
        if(z>eta(i,j)+epsi)
        wval[count] = 0.0;
        
        ++count;
        }
}


