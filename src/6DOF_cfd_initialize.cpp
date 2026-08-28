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

void sixdof_cfd::initialize(lexer *p, fdm *a, ghostcell *pgc)
{
    for (int nb = 0; nb < number6DOF; nb++)
    fb_obj[nb]->initialize_cfd(p, a, pgc);

    if(pchrono)
    pchrono->initialize(p, fb_obj);

    if(number6DOF>1)
    {
        setup(p,a,pgc);

        ALOOP
        fb_union(i,j,k) = 1.0e8;

        for (int nb = 0; nb < number6DOF; nb++)
        {
            fb_obj[nb]->quat_matrices(p);
            fb_obj[nb]->update_position_3D(p,a,pgc,false);
            fb_obj[nb]->add_solid_heaviside(p,a,pgc);
            fb_obj[nb]->union_fb(p,a,fb_union);
        }

        union_solid_field(p,a,pgc);

        pgc->start1(p,a->fbh1,10);
        pgc->start2(p,a->fbh2,11);
        pgc->start3(p,a->fbh3,12);
        pgc->start4(p,a->fbh4,40);

        pgc->gcdf_update(p,a);
    }
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