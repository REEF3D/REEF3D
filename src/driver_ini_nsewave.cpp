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
--------------------------------------------------------------------*/

#include"driver.h"
#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"

void driver::driver_ini_nsewave()
{
    int istart,iend,jstart,jend;
    p->phimean = 0.0;
    p->phiout = 0.0; 
    
    // depth
    if(p->F60>-1.0e20)
    {
    p->phimean=p->F60;
    p->phiout=p->F60;
    p->wd=p->F60;
    }

    // eta plain
    SLICELOOP4
    a->eta(i,j)=0.0;

    // eta slope
    if(p->A251==1)
    SLICELOOP4
    {
    a->eta(i,j)= -p->A251_val*(p->XP[IP]-p->global_xmin);
    }

    // eta box area
    for(int qn=0;qn<p->F72;++qn)
    {
		istart = p->posc_i(p->F72_xs[qn]);
        iend = p->posc_i(p->F72_xe[qn]);

        jstart = p->posc_j(p->F72_ys[qn]);
        jend = p->posc_j(p->F72_ye[qn]);

        SLICELOOP4
        if(i>=istart && i<iend && j>=jstart && j<jend)
        a->eta(i,j) = p->F72_h[qn]-p->wd;
	}

    
    int gcval_phi;
    
    if(p->F50==1)
	gcval_phi=51;

	if(p->F50==2)
	gcval_phi=52;

	if(p->F50==3)
	gcval_phi=53;

	if(p->F50==4)
	gcval_phi=54;
    
    LOOP
    a->phi(i,j,k) = a->eta(i,j) + p->phimean - p->pos_z();
    
    pgc->start4(p,a->phi,gcval_phi);

}