/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
Author: Fabian Knoblauch
--------------------------------------------------------------------*/
#include"VOF_PLIC.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"lexer.h"
#include"gradient.h"

void VOF_PLIC::RKcalcL
(
    fdm *a,
    lexer *p,
    ghostcell* pgc,
    field& uvel,
    field& vvel,
    field& wvel
)
{
    pgc->start4(p,a->vof,1);
    pgc->start1(p,uvel,10);
    pgc->start2(p,vvel,11);
    pgc->start3(p,wvel,12);
    updatePhasemarkers(p,a,pgc);
    starttime=pgc->timer();
    if(p->j_dir>0)
    {
        if(sSweep<5)
            sSweep++;
        else
            sSweep=0;
    }
    else
    {
        if(sSweep<1)
            sSweep++;
        else
            sSweep=0;
    }
    symmetric_scheme2D_FCRK3(a,p,pgc,uvel,vvel,wvel);
    
}