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

#include"pressure_reference.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void pressure_reference::reference_ini(lexer*p, fdm* a, ghostcell *pgc)
{
    double xmin,ymin,zmax;
    int gcinglobal=0;
    
    xmin=ymin=-1.0e10;
    
    gcinglobal=pgc->globalisum(p->gcin_count);
    
    // ini gage location
    if(p->B32==0 && p->B30==1)
    {
    //find active smallest xy inlet location
    
    
    if(gcinglobal>0)
    for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];
    
    xmin = MIN(xmin,p->XP[IP]);
    ymin = MIN(ymin,p->YP[JP]);
    zmax = MAX(zmax,p->ZP[KP]);
    }
    
    if(gcinglobal==0)
    LOOP
    {    
    xmin = MIN(xmin,p->XP[IP]);
    ymin = MIN(ymin,p->YP[JP]);
    zmax = MAX(zmax,p->ZP[KP]);
    }
    
    p->B32_x = pgc->globalmin(xmin);
    p->B32_y = pgc->globalmin(ymin);
    p->B32_z = pgc->globalmax(zmax);
    }
    
}



