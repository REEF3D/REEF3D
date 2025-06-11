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
#include"lexer.h"
#include"ghostcell.h"
void VOF_PLIC::updatePlaneData(lexer* p,fdm* a,ghostcell* pgc ,field& voffield)
{
    pgc->start4(p,voffield,1);
    LOOP
    {
        if(voffield(i,j,k)>=a_thres && voffield(i,j,k)<=w_thres)
        {
            reconstructPlane_alt(a,p,voffield);
            a->nX(i,j,k)=nx(i,j,k);
            a->nY(i,j,k)=ny(i,j,k);
            a->nZ(i,j,k)=nz(i,j,k);
            a->Alpha(i,j,k)=alpha(i,j,k);
        }
        else
        {
            a->nX(i,j,k)=1E06;
            a->nY(i,j,k)=1E06;
            a->nZ(i,j,k)=1E06;
            a->Alpha(i,j,k)=1E06;
        }
    }
    pgc->start4(p,a->nX,1);
    pgc->start4(p,a->nY,1);
    pgc->start4(p,a->nZ,1);
    pgc->start4(p,a->Alpha,1);
}