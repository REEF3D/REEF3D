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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"nhflow_sigma.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"fnpf_ddx_cds2.h"
#include"fnpf_ddx_cds4.h"
#include"fnpf_cds2.h"
#include"fnpf_cds4.h"

#define WLVL (fabs(d->WL(i,j))>1.0e-20?d->WL(i,j):1.0-20)

nhflow_sigma::nhflow_sigma(lexer *p) : nhflow_gradient(p)
{
}

nhflow_sigma::~nhflow_sigma()
{
}

void nhflow_sigma::sigma_coord_ini(lexer *p)
{
    double L, ZN0temp;
    
    L = p->ZN[p->knoz+marge] - p->ZN[0+marge];
    
    ZN0temp = p->ZN[0+marge];
    
    for(k=-marge;k<p->knoz+marge;++k)
    {
    p->ZN[KP] = (p->ZN[KP]-ZN0temp)/L;
    }
}



