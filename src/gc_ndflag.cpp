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
#include"ghostcell.h"
#include"lexer.h"

void ghostcell::ndflag_update(lexer *p)
{
    
    // NDBASEFLAG
    for(i=0;i<p->imax*p->jmax*p->kmax; ++i)
	p->ndbaseflag[i]=p->flag4[i];

	BASELOOP
	{
	    if(p->ndbaseflag[Im1JK]<0)
	    p->ndbaseflag[Im1JK]=9;

	    if(p->ndbaseflag[IJm1K]<0)
	    p->ndbaseflag[IJm1K]=9;

	    if(p->ndbaseflag[IJKm1]<0)
	    p->ndbaseflag[IJKm1]=9;
	}

    BASELOOP
	{

	    if(p->ndbaseflag[Im1JK]>9)
	    if(p->ndbaseflag[IJm1K]>9)
	    p->ndbaseflag[Im1Jm1K]=11;

	    if(p->ndbaseflag[Im1JK]>9)
	    if(p->ndbaseflag[IJKm1]>9)
	    p->ndbaseflag[Im1JKm1]=11;

	    if(p->ndbaseflag[IJm1K]>9)
	    if(p->ndbaseflag[IJKm1]>9)
	    p->ndbaseflag[IJm1Km1]=11;


	    if(p->ndbaseflag[Im1JK]>9)
	    if(p->ndbaseflag[IJm1K]>9)
	    if(p->ndbaseflag[IJKm1]>9)
	    p->ndbaseflag[Im1Jm1Km1]=11;
	}
    
    BASELOOP
    {
    p->ndbaseflag[Im1JK]=11;  
    p->ndbaseflag[IJm1K]=11; 
    p->ndbaseflag[IJKm1]=11;
    
    p->ndbaseflag[Im1Jm1K]=11;
    p->ndbaseflag[Im1JKm1]=11;
    p->ndbaseflag[IJm1Km1]=11;
    
    p->ndbaseflag[Im1Jm1Km1]=11;
    }
}

