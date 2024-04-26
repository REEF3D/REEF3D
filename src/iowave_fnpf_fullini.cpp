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
#include"fdm.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"

void iowave::full_initialize_fnpf(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(p->mpirank==0)
    cout<<"full NWT initialize"<<endl;
    
    // eta
	SLICELOOP4
    if(p->wet[IJ]==1)
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p);

		c->eta(i,j) = wave_eta(p,pgc,xg,yg);

    }
    
    // Fifsf
    SLICELOOP4
    if(p->wet[IJ]==1)
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p);
        
        z = c->eta(i,j);

		c->Fifsf(i,j) = wave_fi(p,pgc,xg,yg,z);
    }

    
    // Fi
    FLOOP
    if(p->wet[IJ]==1)
    {
        xg = xgen(p);
        yg = ygen(p);
        dg = distgen(p);
		db = distbeach(p);
        
        z=p->ZSN[FIJK]-p->phimean;
        
        z = p->ZN[KP]*(c->eta(i,j) + p->wd - c->bed(i,j)) + c->bed(i,j)-p->phimean;
        
        c->Fi[FIJK] = wave_fi(p,pgc,xg,yg,z);
      
    }
    
    SLICELOOP4
    c->WL(i,j) = c->eta(i,j) + p->wd - c->bed(i,j);

}

