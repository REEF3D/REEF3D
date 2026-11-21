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

#include"iowave.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"

void iowave::active_beach2D(lexer *p, fdm2D* b, ghostcell* pgc, slice &P, slice &Q, slice &bed, slice &eta)
{
    int aa,bb;

    // P at the beach faces (B 99 3: shallow-water radiation velocity sqrt(g/d) eta; B 99 4: P carried out of the
    // domain with the wave speed sqrt(g h)), three ghost faces
    for(n=0;n<p->gcslawa1_count;++n)
    {
    i=p->gcslawa1[n][0];
    j=p->gcslawa1[n][1];

    const int side = p->gcslawa1[n][2];
    const double eta_R = eta(i+(side==4 ? 1 : 0),j);

    aa = (side==1) ? -1 : (side==4) ? 1 : 0;
    bb = (side==3) ? -1 : (side==2) ? 1 : 0;

        if(p->B99==3)
        {
        const double Uc = eta_R*sqrt(9.81/p->wd);

        P(i+1*aa,j+1*bb) = Uc;
        P(i+2*aa,j+2*bb) = Uc;
        P(i+3*aa,j+3*bb) = Uc;
        }

        if(p->B99==4)
        {
        const double dfx4 = (P(i,j)-P(i-1,j))/p->DXM;

        for(int q=1;q<=3;++q)
        P(i+q*aa,j+q*bb) = P(i,j) - p->dt*sqrt(9.81*b->hp(i,j))*dfx4;
        }
    }

    // Q at the beach faces: no velocity in the ghost faces
    for(n=0;n<p->gcslawa2_count;++n)
    {
    i=p->gcslawa2[n][0];
    j=p->gcslawa2[n][1];

    const int side = p->gcslawa2[n][2];

    aa = (side==1) ? -1 : (side==4) ? 1 : 0;
    bb = (side==3) ? -1 : (side==2) ? 1 : 0;

    b->Q(i+1*aa,j+1*bb) = 0.0;
    b->Q(i+2*aa,j+2*bb) = 0.0;
    b->Q(i+3*aa,j+3*bb) = 0.0;
    }
}
