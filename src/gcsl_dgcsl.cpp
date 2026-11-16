/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
#include"slice.h"

void ghostcell::dgcslpol1(lexer* p,slice& f)
{

    for(n=0;n<p->dgcsl1_count;++n)
    {
    i=p->dgcsl1[n][0];
    j=p->dgcsl1[n][1];


    if(p->dgcsl1[n][2]==1)
    f(i-1,j-1) = f(i,j);
    
    if(p->dgcsl1[n][2]==2)
    f(i+1,j-1) = f(i,j);
    
    if(p->dgcsl1[n][2]==3)
    f(i+1,j+1) = f(i,j);
    
    if(p->dgcsl1[n][2]==4)
    f(i-1,j+1) = f(i,j);
    }
}

void ghostcell::dgcslpol2(lexer* p,slice& f)
{

    for(n=0;n<p->dgcsl2_count;++n)
    {
    i=p->dgcsl2[n][0];
    j=p->dgcsl2[n][1];


    if(p->dgcsl2[n][2]==1)
    f(i-1,j-1) = f(i,j);
    
    if(p->dgcsl2[n][2]==2)
    f(i+1,j-1) = f(i,j);
    
    if(p->dgcsl2[n][2]==3)
    f(i+1,j+1) = f(i,j);
    
    if(p->dgcsl2[n][2]==4)
    f(i-1,j+1) = f(i,j);
    }
}

void ghostcell::dgcslpol4(lexer* p,slice& f)
{

    for(n=0;n<p->dgcsl4_count;++n)
    {
    i=p->dgcsl4[n][0];
    j=p->dgcsl4[n][1];


    if(p->dgcsl4[n][2]==1)
    f(i-1,j-1) = f(i,j);
    
    if(p->dgcsl4[n][2]==2)
    f(i+1,j-1) = f(i,j);
    
    if(p->dgcsl4[n][2]==3)
    f(i+1,j+1) = f(i,j);
    
    if(p->dgcsl4[n][2]==4)
    f(i-1,j+1) = f(i,j);
    }
}
