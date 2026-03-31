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

#include"vrans_nhflow_f.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

vrans_nhflow_f::vrans_nhflow_f(lexer *p, fdm_nhf *d, ghostcell *pgc) : nhflow_geometry(p,d,pgc), Cval(p->B264)
{
	p->Darray(NPOR,p->imax*p->jmax*(p->kmax+2));
    p->Darray(DPOR,p->imax*p->jmax*(p->kmax+2));
    p->Darray(APOR,p->imax*p->jmax*(p->kmax+2));
    p->Darray(BPOR,p->imax*p->jmax*(p->kmax+2));
}

vrans_nhflow_f::~vrans_nhflow_f()
{
}

void vrans_nhflow_f::update(lexer *p, fdm_nhf *d, ghostcell *pgc, int val)
{
    ray_cast(p, d, pgc, d->PORSTRUC);
    reini_RK2(p, d, pgc, d->PORSTRUC);
    

    LOOP
    {
    H = Hporface(p,d,0,0,0);   
    
    d->POR[IJK]     = H*p->B201_n;
	d->PORPART[IJK] = H*p->B201_d50;
	APOR[IJK]       = H*p->B201_alpha;
	BPOR[IJK]       = H*p->B201_beta;
    }
    
    pgc->start4V(p,d->POR,1);
    pgc->start4V(p,d->PORPART,1);
    pgc->start4V(p,APOR,1);
    pgc->start4V(p,BPOR,1);
    
}

