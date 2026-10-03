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

#include"bc_noflux.h"
#include"lexer.h"

void bc_noflux_mask(lexer *p, std::vector<int> &mask, int what)
{
    mask.assign(size_t(p->imax)*size_t(p->jmax)*size_t(p->kmax),0);

    const int fixed_scalar[7] = {0,p->H61,p->H62,p->H63,p->H64,p->H65,p->H66};

    for(int n=0;n<p->gcb4_count;++n)
    {
        int bc = p->gcb4[n][4];
        int cs = p->gcb4[n][3];

        if(cs<1 || cs>6)
        continue;

        if(bc_periodic_face(p,p->gcb4[n][0],p->gcb4[n][1],p->gcb4[n][2],cs))
        continue;

        bool noflux=false;

        if((what & BC_NOFLUX_WALLS) && (bc==3 || bc==5 || bc==21 || bc==22))
        noflux=true;

        if((what & BC_NOFLUX_INFLOW) && (bc==1 || bc==6))
        noflux=true;

        // patchBC codes: hundreds = velocity, tens = pressure (1 Neumann, 2 Dirichlet)
        if((what & BC_NOFLUX_INFLOW) && bc>=100 && bc<=222 && bc/100==2 && (bc/10)%10==1)
        noflux=true;

        if((what & BC_NOFLUX_SCALAR) && fixed_scalar[cs]==1)
        noflux=false;

        if(noflux)
        {
        int i=p->gcb4[n][0];
        int j=p->gcb4[n][1];
        int k=p->gcb4[n][2];

        mask[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin] |= (1<<(cs-1));
        }
    }
}

bool bc_periodic_face(lexer *p, int i, int j, int k, int cs)
{
    if(cs==1 && p->periodic1==1 && i+p->origin_i==0)
    return true;

    if(cs==4 && p->periodic1==1 && i+p->origin_i==p->gknox-1)
    return true;

    if(cs==3 && p->periodic2==1 && j+p->origin_j==0)
    return true;

    if(cs==2 && p->periodic2==1 && j+p->origin_j==p->gknoy-1)
    return true;

    if(cs==5 && p->periodic3==1 && k+p->origin_k==0)
    return true;

    if(cs==6 && p->periodic3==1 && k+p->origin_k==p->gknoz-1)
    return true;

    return false;
}
