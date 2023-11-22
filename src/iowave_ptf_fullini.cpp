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
/*

#include"iowave.h"
#include"lexer.h"
#include"fdm_ptf.h"
#include"ghostcell.h"

void iowave::full_initialize_ptf(lexer *p, fdm_ptf *e, ghostcell *pgc)
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

		e->eta(i,j) = wave_eta(p,pgc,xg,yg);

    }
    
    // Fifsf
    SLICELOOP4
    if(p->wet[IJ]==1)
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p);
        
        z = e->eta(i,j);

		e->Fifsf(i,j) = wave_fi(p,pgc,xg,yg,z);
    }

    
    // Fi
    FLUIDLOOP
    if(p->wet[IJ]==1)
    {
        xg = xgen(p);
        yg = ygen(p);
        dg = distgen(p);
		db = distbeach(p);
        
        z=p->ZN[KP]-p->phimean;
        
        z = p->ZN[KP]*(e->eta(i,j) + p->wd - e->bed(i,j)) + e->bed(i,j)-p->phimean;
        
        e->Fi[IJK] = wave_fi(p,pgc,xg,yg,z);
      
    }
    
    SLICELOOP4
    e->WL(i,j) = e->eta(i,j) + p->wd - e->bed(i,j);

}

*/