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
#include"patchBC_codes.h"

// ghost-cell label of u (kernel chosen in gcdistro1)
// gcv: 1/10 u, 7 projection, 14 / 110 / 114 / 117 special treatments, 50 Neumann; cs: side of the face
int ghostcell::gceval1(lexer *p, int gcv, int bc, int cs)
{
    const bool para = (cs==2||cs==3||cs==5||cs==6);          // face parallel to u
    const bool orth = (cs==1||cs==4);                        // face normal to u
    const bool wall = (bc==21||bc==22||(bc==7&&awa_lable==0)); // walls; beach without active absorption
    const bool walltopo = (bc==21||bc==22||bc==5);
    const bool patch_out = patch_outlet(bc);                 // patch outlet: zero gradient, inlet: set by patchBC
    const bool outflow = (bc==2 && gclabel_outflow==1);

    switch(gcv)
    {
    case 50:
        return 4;

    case 1:
    case 10:
        if(para)
        {
            if(wall)                return gclabel_u;
            if(bc==5)               return gclabel_utopo;
            if(outflow || bc==3)    return 4;
            if(bc==9 && cs==6)      return 4;
        }
        if(orth)
        {
            if(walltopo || wall || bc==3) return gclabel_u_orth;
            if(bc==6)               return gclabel_u_in;
            if(outflow)             return gclabel_u_out;
            if((bc==7||bc==8) && gclabel_outflow==1 && p->I10==1) return 4;
        }
        if(patch_out)               return 4;
        return 0;

    case 7:
        if(orth && (walltopo || wall)) return gclabel_vel;
        if(orth && bc==6)           return gclabel_u_in;
        if(patch_out)               return 4;
        return 0;

    case 14:
        return walltopo ? 4 : 0;

    case 17:
        if(bc==3 && para)           return 4;
        if(bc==3 && orth)           return gclabel_u_orth;
        if(bc==9 && cs==6)          return 4;
        return 0;

    case 110:
        return (para && (wall || bc==5)) ? 5 : 0;

    case 114:
        if(para && wall)            return gclabel_u;
        if(para && bc==5)           return gclabel_utopo;
        return 0;

    case 117:
        if(para && wall)            return 4;
        if(para && bc==5)           return gclabel_utopo;
        return 0;
    }

    return 0;
}

void ghostcell::gcdistro1(lexer *p,field& f, int ii, int jj, int kk, int nn, double dist,  int gcv, int bc, int cs)
{
    i=ii;
	j=jj;
	k=kk;
	n=nn;

	bc_label=gceval1(p,gcv,bc,cs);

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

