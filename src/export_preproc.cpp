/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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

#include"export.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void exportfile::preproc(lexer *p, fdm *a, ghostcell *pgc)
{
    xs = p->originx;
    ys = p->originy;
    zs = p->originz;
    
    xe = (double(p->knox)+0.5)*p->DXM + p->originx;
    ye = (double(p->knoy)+0.5)*p->DXM + p->originy;
    ze = (double(p->knoz)+0.5)*p->DXM + p->originz;
    /*
    LOOP
    eta[i][j]=p->global_zmin-1.0e20;
    
    // eta 
    LOOP
    if(a->phi(i,j,k)>=0.0 && a->phi(i,j,k+1)<0.0)
    eta[i][j] = MAX(eta[i][j], -(a->phi(i,j,k)*p->DXM)/(a->phi(i,j,k+1)-a->phi(i,j,k)) + p->pos_z());
    
    pgc->verticalmax(p,a,eta);*/
    
    
}
