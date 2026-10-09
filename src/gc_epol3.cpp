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

// ghost-cell label of w (kernel chosen in gcdistro3)
// gcv: 3/12 w, 9 projection, 16 / 112 / 116 / 119 special treatments, 17 omega_sig, 19, 50 Neumann;
// cs: side of the face. The bottom face (cs 5) of walls is treated only for A10 = 5 and 6;
// the lid (3) is a wall for CFD, Neumann for FNPF (A10 3), not set for NHFLOW (A10 5)
int ghostcell::gceval3(lexer *p, int gcv, int bc, int cs)
{
    const bool para = (cs==1||cs==2||cs==3||cs==4);          // face parallel to w
    const bool orth = (cs==5||cs==6);                        // face normal to w
    const bool wall = (bc==21||bc==22||(bc==7&&awa_lable==0));
    const bool walltopo = (bc==21||bc==22||bc==5);
    const bool patch = (bc==111||bc==112||bc==121||bc==122);
    const bool outflow = (bc==2 && gclabel_outflow==1);

    switch(gcv)
    {
    case 50:
        return 4;

    case 3:
    case 12:
        if(para)
        {
            if(wall)                return gclabel_w;
            if(bc==5)               return gclabel_wtopo;
            if(outflow || bc==3)    return 4;
        }
        if(orth)
        {
            if((walltopo || wall) && (cs==6 || p->A10==5 || p->A10==6)) return gclabel_w_orth;
            if(outflow)             return gclabel_w_out;
            if(bc==3)               return p->A10==3 ? 4 : (p->A10==5 ? 0 : gclabel_w_orth);
            if(bc==9 && cs==6)      return 4;
        }
        if(bc==6)                   return gclabel_w_in;
        if(patch)                   return 4;
        return 0;

    case 9:
        if(orth && (walltopo || wall)) return gclabel_vel;
        if(bc==6)                   return gclabel_w_in;
        if(patch)                   return 4;
        return 0;

    case 16:
        return walltopo ? 4 : 0;

    case 17:
        if(bc==21 && cs==5)         return 5;
        if(bc==3 && cs==6)          return 0;
        return 4;

    case 19:
        if(bc==3 && para)           return 4;
        if(bc==3 && orth)           return p->A10==3 ? 4 : (p->A10==5 ? 0 : gclabel_w_orth);
        if(bc==9 && cs==6)          return 4;
        return 0;

    case 112:
        return (para && (wall || bc==5)) ? 5 : 0;

    case 116:
        if(para && wall)            return gclabel_w;
        if(para && bc==5)           return gclabel_wtopo;
        return 0;

    case 119:
        if(para && wall)            return 4;
        if(para && bc==5)           return gclabel_wtopo;
        return 0;
    }

    return 0;
}

void ghostcell::gcdistro3(lexer *p,field& f, int ii, int jj, int kk, int nn, double dist,  int gcv, int bc, int cs)
{
    i=ii;
	j=jj;
	k=kk;
	n=nn;

	bc_label=gceval3(p,gcv,bc,cs);

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

