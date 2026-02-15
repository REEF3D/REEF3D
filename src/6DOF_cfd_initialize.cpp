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
Authors: Tobias Martin, Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_cfd.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<algorithm>

void sixdof_cfd::initialize(lexer *p, fdm *a, ghostcell *pgc)
{
    for (int nb = 0; nb < number6DOF; nb++)
    {
    fb_obj[nb]->initialize_cfd(p, a, pgc);
    
    fb_combine(p,a,nb);
    }
    
    fb_finish(p,a,pgc);
}

void sixdof_cfd::initialize(lexer *p, fdm2D *b, ghostcell *pgc)
{
}

void sixdof_cfd::initialize(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
}

void sixdof_cfd::ini(lexer *p, ghostcell *pgc)
{
}

void sixdof_cfd::setup(lexer *p, fdm *a, ghostcell *pgc)
{
    // the weights of the solids and the topography (ghostcell::solid_forcing of this stage)
    if(sf[0]!=nullptr)
    {
    const int n = p->imax*p->jmax*p->kmax;
    std::copy(a->fbh1.data(), a->fbh1.data()+n, sf[0]->data());
    std::copy(a->fbh2.data(), a->fbh2.data()+n, sf[1]->data());
    std::copy(a->fbh3.data(), a->fbh3.data()+n, sf[2]->data());
    }
    
    // Reset heaviside field
    ULOOP
    a->fbh1(i,j,k) = 0.0;

    VLOOP
    a->fbh2(i,j,k) = 0.0;
    
    WLOOP
    a->fbh3(i,j,k) = 0.0;

    LOOP
    a->fbh4(i,j,k) = 0.0;

    pgc->start1(p,a->fbh1,10);
    pgc->start2(p,a->fbh2,11);
    pgc->start3(p,a->fbh3,12);
    pgc->start4(p,a->fbh4,40);
}