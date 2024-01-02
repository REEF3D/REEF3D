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
#include"fdm2D.h"
#include"ghostcell.h"

void iowave::full_initialize2D(lexer *p, fdm2D *b, ghostcell *pgc)
{
	 double uval,vval;
    double deltaz;
    
    SLICELOOP1
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p);
        

        
        deltaz = (0.5*(b->eta(i,j)+b->eta(i+1,j)) + p->wd - 0.5*(b->bed(i,j)+b->bed(i+1,j)))/(double(p->B160));
        
        uval=0.0;
        z=-p->wd;
        for(int qn=0;qn<p->B160;++qn)
        {
        uval += wave_u(p,pgc,xg,yg,z);
        
        z+=deltaz;
        }
        uval/=double(p->B160);
        

		b->P(i,j) = uval;
	}
	
    
    SLICELOOP2
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p);


        deltaz = (0.5*(b->eta(i,j)+b->eta(i,j+1)) + p->wd - 0.5*(b->bed(i,j)+b->bed(i,j+1)))/(double(p->B160));
        
        vval=0.0;
        z=-p->wd;
        for(int qn=0;qn<p->B160;++qn)
        {
        vval += wave_v(p,pgc,xg,yg,z);
        
        z+=deltaz;
        }
        vval/=double(p->B160);
        
        b->Q(i,j) = vval;
	}
	
	// eta
	SLICELOOP4
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p);

		b->eta(i,j) = wave_eta(p,pgc,xg,yg);

    }
}
