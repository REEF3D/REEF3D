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

#include"nhflow_forcing.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"nhflow_reinidisc_fsf.h"

nhflow_forcing::nhflow_forcing(lexer *p) : epsi(1.6)
{
    p->Iarray(IO,p->imax*p->jmax*(p->kmax+2));
    p->Iarray(CL,p->imax*p->jmax*(p->kmax+2));
    p->Iarray(CR,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(FX,p->imax*p->jmax*(p->kmax+2));
    p->Darray(FY,p->imax*p->jmax*(p->kmax+2));
    p->Darray(FZ,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(FRK1,p->imax*p->jmax*(p->kmax+2));
    p->Darray(dt,p->imax*p->jmax*(p->kmax+2));
    p->Darray(L,p->imax*p->jmax*(p->kmax+2));
    
    prdisc = new nhflow_reinidisc_fsf(p);
}

nhflow_forcing::~nhflow_forcing()
{
}

void nhflow_forcing::forcing(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha, double *U, double *V, double *W)
{
    // update direct forcing function
    ray_cast(p, d, pgc);
    reini_RK2(p, d, pgc, d->SOLID);
    
    // update Heaviside
    
    // add forcing term to RHS
}

void nhflow_forcing::forcing_ini(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    LOOP
    p->ZSP[IJK]  = p->ZP[KP]*d->WL(i,j) + d->bed(i,j);
    
    objects_create(p, pgc);
    
    geometry_refinement(p,pgc);
    
    ray_cast(p, d, pgc);
    
    reini_RK2(p, d, pgc, d->SOLID);
}