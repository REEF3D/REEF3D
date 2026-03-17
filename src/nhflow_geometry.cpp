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

#include"nhflow_geometry.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"6DOF.h"
#include"nhflow_reinidisc_fsf.h"

nhflow_geometry::nhflow_geometry(lexer *p, fdm_nhf *d, ghostcell *pgc) : epsi(1.6)
{
    forcing_flag=0;
    solid_flag=0;
    floating_flag=0;
    dlm_flag=0;
    
    if(p->A581>0 || p->A583>0 || p->A584>0   || p->A585>0  || p->A586>0 || p->A587>0 || p->A588>0 || p->A589>0 || p->A590>0)
    {   
    forcing_flag=1;
    solid_flag=1;
    }
    
    if(p->X10>0)
    {
    forcing_flag=1;
    floating_flag=1;
    }
    
    if(p->A599==1)
    {
    dlm_flag=1;
    forcing_flag=0;
    solid_flag=0;
    floating_flag=0;
    }
    
    // ----
    if(forcing_flag==1)
    {
    p->Iarray(IO,p->imax*p->jmax*(p->kmax+2));
    p->Iarray(CL,p->imax*p->jmax*(p->kmax+2));
    p->Iarray(CR,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(FRK1,p->imax*p->jmax*(p->kmax+2));
    p->Darray(dt,p->imax*p->jmax*(p->kmax+2));
    p->Darray(L,p->imax*p->jmax*(p->kmax+2));

    prdisc = new nhflow_reinidisc_fsf(p);
    }
}

nhflow_geometry::~nhflow_geometry()
{
}

