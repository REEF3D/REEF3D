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
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include<math.h>

// ghost-cell label of the topo / solid / fb / level-set fields in start4a (kernel chosen in gcdistro4a):
// 4 neumann, 3 extend, 0 nothing; cs: side of the face (sign removed)
int ghostcell::gceval4a(lexer *p, int gcv, int bc, int cs)
{
    switch(gcv)
    {
    // topo, topo for bedload
    case 151:
    case 152:
    case 153:
    case 161:
    case 162:
    case 163:
        if(bc==21||bc==22||bc==5||bc==3||bc==6||bc==7||bc==8) return 4;
        if(bc==2 && (gcv==151||gcv==161))           return 4;
        if(bc==1 && (gcv==152||gcv==162))           return 4;
        return 0;

    case 150:
    case 154:
        return 4;

    case 159:
        return 3;

    // floating body
    case 50:
        return 4;

    // level set
    case 51:
    case 52:
    case 53:
    case 54:
        if(bc==21||bc==22||bc==5||bc==41||bc==3||bc==6||bc==7||bc==8||bc==9) return 4;
        if(bc==1 && (gcv==52||gcv==54))             return 4;
        if(bc==2 && (gcv==51||gcv==54))             return 4;
        return 0;

    // porosity
    case 1:
        return 4;
    }

    return 0;
}

void ghostcell::gcdistro4a(lexer *p,field& f, int ii, int jj, int kk, int nn, double dist,  int gcv, int bc, int cs)
{
    i=ii;
	j=jj;
	k=kk;
	n=nn;
	
	cs = fabs(cs);
    

	bc_label=gceval4a(p,gcv,bc,cs);

    switch(bc_label)
    {
    case 4: neumann(f,gcv,bc,cs); break;
    case 3: extend(p,f,dist,gcv,bc,cs); break;
    }
}


