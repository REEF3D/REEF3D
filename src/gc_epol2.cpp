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

// ghost-cell label of v (kernel chosen in gcdistro2)
// gcv: 2/11 v, 8 projection, 15 / 111 / 115 / 118 special treatments, 18, 50 Neumann; cs: side of the face
// differences to u and w: wavegen (6) and beach (7, also with active absorption) are walls for the
// parallel faces, the beach also for the orthogonal faces; the lid (3) uses label 1 instead of gclabel_v_orth
int ghostcell::gceval2(lexer *p, int gcv, int bc, int cs)
{
    const bool para = (cs==1||cs==4||cs==5||cs==6);          // face parallel to v
    const bool orth = (cs==2||cs==3);                        // face normal to v
    const bool wall = (bc==21||bc==22||bc==7||bc==6);
    const bool walltopo = (bc==21||bc==22||bc==5);
    const bool patch = (bc==111||bc==112||bc==121||bc==122);
    const bool outflow = (bc==2 && gclabel_outflow==1);

    switch(gcv)
    {
    case 50:
        return 4;

    case 2:
    case 11:
        if(para)
        {
            if(wall)                return gclabel_v;
            if(bc==5)               return gclabel_vtopo;
            if(outflow || bc==3)    return 4;
            if(bc==9 && cs==6)      return 4;
        }
        if(orth)
        {
            if(walltopo || bc==7)   return gclabel_v_orth;
            if(outflow)             return gclabel_v_out;
            if(bc==3)               return 1;
        }
        if(bc==6)                   return gclabel_v_in;
        if(patch)                   return 4;
        return 0;

    case 8:
        if(orth && (walltopo || bc==7)) return gclabel_vel;
        if(bc==6)                   return gclabel_v_in;
        if(patch)                   return 4;
        return 0;

    case 15:
        return walltopo ? 4 : 0;

    case 18:
        if(bc==3 && para)           return 4;
        if(bc==3 && orth)           return 1;
        if(bc==9 && cs==6)          return 4;
        return 0;

    case 111:
        return (para && (wall || bc==5)) ? 5 : 0;

    case 115:
        if(para && wall)            return gclabel_v;
        if(para && bc==5)           return gclabel_vtopo;
        return 0;

    case 118:
        return (para && (wall || bc==5)) ? 4 : 0;
    }

    return 0;
}


void ghostcell::gcdistro2(lexer *p,field& f, int ii, int jj, int kk, int nn, double dist,  int gcv, int bc, int cs)
{
    i=ii;
	j=jj;
	k=kk;
	n=nn;

	bc_label=gceval2(p,gcv,bc,cs);

    switch(bc_label)
    {
    case 1: dirichlet_ortho(p,f,dist,gcv,bc,cs); break;
    case 2: dirichlet_para(p,f,dist,gcv,bc,cs); break;
    case 4: neumann(f,gcv,bc,cs); break;
    case 5: noslip(f,dist,gcv,bc,cs); break;
    case 6: outflow(p,f,gcv,bc,cs); break;
    case 11: dirichlet_ortho_reflect(p,f,dist,gcv,bc,cs); break;
    case 12: dirichlet_para_reflect(p,f,dist,gcv,bc,cs); break;
    }
}
