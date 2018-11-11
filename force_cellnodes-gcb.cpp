/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"force.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void force::cellnodes_gcb(lexer* p, fdm *a, ghostcell *pgc)
{
    const double eps = 1.0e-5*p->dx;
    int cn[4][4];
    int r,s;
    int count, fcount;

    fcount=0;

    // identify active and unactive cell nodes
    GC4LOOP
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];

		if(i>=is && i<=ie && j>=js && j<=je && k>=ks && k<=ke && p->gcb4[n][4]==21)
		{		   

			fid[fcount][0] = i;
			fid[fcount][1] = j;
			fid[fcount][2] = k;
			
			fx[fcount] = p->XP[IP];
			fy[fcount] = p->YP[JP];
			fz[fcount] = p->ZP[KP];

			
			if(p->gcb4[n][3]==1)
			{
			fn[fcount][0] = 1.0;
			fn[fcount][1] = 0.0;
			fn[fcount][2] = 0.0;
			fx[fcount] -= 0.75*p->DXP[IP];
            
            farea[fcount] = p->DYN[JP]*p->DZP[KP]; 
			}
			
			if(p->gcb4[n][3]==4)
			{
			fn[fcount][0] = -1.0;
			fn[fcount][1] = 0.0;
			fn[fcount][2] = 0.0;
			fx[fcount] += 0.75*p->DXP[IP];
            
            farea[fcount] = p->DYN[JP]*p->DZP[KP];
			}
			
			if(p->gcb4[n][3]==3)
			{
			fn[fcount][0] = 0.0;
			fn[fcount][1] = 1.0;
			fn[fcount][2] = 0.0;
			fy[fcount] -= 0.75*p->DYP[JP];
            
            farea[fcount] = p->DXN[IP]*p->DZP[KP];
			}
			
			if(p->gcb4[n][3]==2)
			{
			fn[fcount][0] = 0.0;
			fn[fcount][1] = -1.0;
			fn[fcount][2] = 0.0;
			fy[fcount] += 0.75*p->DYP[JP];
            
            farea[fcount] = p->DXN[IP]*p->DZP[KP];
			}
			
			if(p->gcb4[n][3]==5)
			{
			fn[fcount][0] = 0.0;
			fn[fcount][1] = 0.0;
			fn[fcount][2] = 1.0;
			fz[fcount] -= 0.75*p->DZP[KP];
            
            farea[fcount] = p->DXN[IP]*p->DYP[JP];
			}
			
			if(p->gcb4[n][3]==6)
			{
			fn[fcount][0] = 0.0;
			fn[fcount][1] = 0.0;
			fn[fcount][2] = -1.0;
			fz[fcount] += 0.75*p->DZP[KP];
            
            farea[fcount] = p->DXN[IP]*p->DYP[JP];
			}
			++fcount;
		}
    }
}


