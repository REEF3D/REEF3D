/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"vrans_nhflow.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

vrans_nhflow::vrans_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc) : nhflow_geometry(p,d,pgc), Cval(p->B264)
{
	p->Darray(NPOR,p->imax*p->jmax*(p->kmax+2));
    p->Darray(DPOR,p->imax*p->jmax*(p->kmax+2));
    p->Darray(APOR,p->imax*p->jmax*(p->kmax+2));
    p->Darray(BPOR,p->imax*p->jmax*(p->kmax+2));
}

vrans_nhflow::~vrans_nhflow()
{
}

