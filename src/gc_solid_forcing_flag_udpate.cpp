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
Authors: Hans Bihs
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"

void ghostcell::solid_forcing_flag_update(lexer *p, fdm *a)
{
    // Update DF
    LOOP
    p->DF[IJK]=1;
    
    if(p->topoforcing>0 && p->solidread>0)
    LOOP
    if(a->solid(i,j,k)<0.0 || a->topo(i,j,j)<0.0)
    p->DF[IJK]=-1;
    
    if(p->topoforcing>0 && p->solidread==0)
    LOOP
    if(a->topo(i,j,j)<0.0)
    p->DF[IJK]=-1;
    
    if(p->topoforcing==0 && p->solidread>0)
    LOOP
    if(a->solid(i,j,k)<0.0)
    p->DF[IJK]=-1;
    
    
    startintV(p,p->DF,1);
}