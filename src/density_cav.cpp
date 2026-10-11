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
Author: Ahmet Soydan
--------------------------------------------------------------------*/


#include"density_cav.h"
#include"cavitation.h"
#include"lexer.h"
#include"fdm.h"

density_cav::density_cav(lexer* p)
{
}

density_cav::~density_cav()
{
}

double density_cav::roface(lexer *p, fdm *a, int aa, int bb, int cc)
{
    if(cavitation::active==nullptr)
    return p->W1;

    const double al = 0.5*(cavitation::active->alpha(i,j,k) + cavitation::active->alpha(i+aa,j+bb,k+cc));

    return cavitation::active->rho_mix(p,al);
}
