/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"ghostcell.h"
#include"lexer.h"
#include"fdm_nhf.h"

void ghostcell::gciobc_update(lexer *p, fdm_nhf *d)
{
    /*
    LOOP
    p->IO[IJK]=0;
    
    LOOP
    if(d->SOLID[IJK]>0 && d->FB[IJK]>0)
    {
    if(p->flag4[Im1JK]<0)
    p->IO[Im1JK]=-10;
    
    if(p->flag4[Ip1JK]<0)
    p->IO[Ip1JK]=-10;
    
    if(p->flag4[IJm1K]<0)
    p->IO[IJm1K]=-10;
    
    if(p->flag4[IJp1K]<0)
    p->IO[IJp1K]=-10;
    
    if(p->flag4[IJKm1]<0)
    p->IO[IJKm1]=-10;
    
    if(p->flag4[IJKp1]<0)
    p->IO[IJKp1]=-10;
    
    
    if(d->SOLID[Im1JK]<0)
    p->IO[Im1JK]=-10;
    
    if(d->SOLID[Ip1JK]<0)
    p->IO[Ip1JK]=-10;
    
    if(d->SOLID[IJm1K]<0)
    p->IO[IJm1K]=-10;
    
    if(d->SOLID[IJp1K]<0)
    p->IO[IJp1K]=-10;
    
    if(d->SOLID[IJKm1]<0)
    p->IO[IJKm1]=-10;
    
    if(d->SOLID[IJKp1]<0)
    p->IO[IJKp1]=-10;
    
    
    
    if(d->FB[Im1JK]<0)
    p->IO[Im1JK]=-10;
    
    if(d->FB[Ip1JK]<0)
    p->IO[Ip1JK]=-10;
    
    if(d->FB[IJm1K]<0)
    p->IO[IJm1K]=-10;
    
    if(d->FB[IJp1K]<0)
    p->IO[IJp1K]=-10;
    
    if(d->FB[IJKm1]<0)
    p->IO[IJKm1]=-10;
    
    if(d->FB[IJKp1]<0)
    p->IO[IJKp1]=-10;
    }
    
    startintV(p,p->IO,1);*/
    
}