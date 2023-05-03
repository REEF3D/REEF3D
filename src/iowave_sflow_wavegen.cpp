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
#include"fdm2D.h"
#include"ghostcell.h"

void iowave::wavegen2D(lexer *p, fdm2D* b, ghostcell* pgc, slice &P, slice &Q, slice &bed, slice &eta)
{
    double uval,vval,wval;
    double deltaz;
    
	for(n=0;n<p->gcslin_count;n++)
    {
		i=p->gcslin[n][0];
		j=p->gcslin[n][1];

        x=xgen(p);
        y=ygen(p);
        x1=xgen1(p);
        y2=ygen2(p);
		
        // u
		deltaz = (0.5*(eta(i,j)+eta(i+1,j)) + p->wd - 0.5*(bed(i,j)+bed(i+1,j)))/(double(p->B160));
        
        uval=0.0;
        z=-p->wd;
        for(int qn=0;qn<=p->B160;++qn)
        {
        uval += wave_u(p,pgc,x1,y,z);
        
        z+=deltaz;
        }
        uval/=double(p->B160+1);
        
        // v
        deltaz = (0.5*(eta(i,j)+eta(i,j+1)) + p->wd - 0.5*(bed(i,j)+bed(i,j+1)))/(double(p->B160));
        
        vval=0.0;
        z=-p->wd;
        for(int qn=0;qn<=p->B160;++qn)
        {
        vval += wave_v(p,pgc,x,y2,z);
        
        z+=deltaz;
        }
        vval/=double(p->B160+1);
        
        // w
        z=eta(i,j);

        wval = wave_w(p,pgc,xg,yg,z);
        

        P(i-1,j)=uval+p->Ui;
        P(i-2,j)=uval+p->Ui;
        P(i-3,j)=uval+p->Ui;
        

        Q(i-1,j)=vval;
        Q(i-2,j)=vval;
        Q(i-3,j)=vval;
        
        b->ws(i-1,j)=wval;
        b->ws(i-2,j)=wval;
        b->ws(i-3,j)=wval;
        
        eta(i,j)  =wave_eta(p,pgc,x,y);
        eta(i-1,j)=wave_eta(p,pgc,x,y);
        eta(i-2,j)=wave_eta(p,pgc,x,y);
        eta(i-3,j)=wave_eta(p,pgc,x,y);
    }	
}
