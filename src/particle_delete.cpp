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

#include"particle_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void particle_f::random_delete(lexer* p, fdm* a, ghostcell* pgc)
{
	double lsc, maxpos;
	double ppcell_coeff=1.0;
	
	int qn;
	
	// POS
	for(qn=0;qn<2;++qn)
    for(n=0;n<posactive;++n)
    if(posflag[n]==1)
    {
		if(qn==0)
		ppcell_coeff = 1.25;
		
		if(qn==1)
		ppcell_coeff = 2.25;
		
        i=int((pos[n][0])/dx);
        j=int((pos[n][1])/dy);
        k=int((pos[n][2])/dz);
		
		lsc = a->phi(i,j,k);
		
		if(lsc<0.5*p->DXM && lsc>-0.5*p->DXM)
		maxpos = (0.5 + lsc/p->DXM)*double(ppcell)*ppcell_coeff;
		
		if(lsc>=0.5*p->DXM)
		maxpos = double(ppcell)*ppcell_coeff;
		
		if(lsc<=-0.5*p->DXM)
		maxpos = 0.0;
		
		if(posnum(i,j,k)>maxpos)		
		if(pos[n][3]>pos[n][4]+0.1*rmin || qn==1)
		{
		++pcount;
		posflag[n]=0;
        posmem[pcount]=n;
        ++removed;
		posnum(i,j,k)-=1.0;
		}    
    }
}
