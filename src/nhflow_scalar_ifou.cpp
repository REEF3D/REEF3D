/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"nhflow_scalar_ifou.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"nhflow_scalar_advec_CDS2.h"

nhflow_scalar_ifou::nhflow_scalar_ifou(lexer *p)
{

    padvec = new nhflow_scalar_advec_CDS2(p);

}

nhflow_scalar_ifou::~nhflow_scalar_ifou()
{
}

void nhflow_scalar_ifou::start(lexer* p, fdm_nhf *d, double *F, int ipol, double *U, double *V, double *W)
{
    count=0;
    LOOP
    {
    udir=vdir=wdir=0.0;
    
    padvec->uadvec(ipol,U,ivel1,ivel2);
    padvec->vadvec(ipol,V,jvel1,jvel2);
    padvec->wadvec(ipol,W,kvel1,kvel2);

	if(0.5*(ivel1+ivel2)>=0.0)
    udir=1.0;
    
    if(0.5*(jvel1+jvel2)>=0.0)
    vdir=1.0;
    
    if(0.5*(kvel1+kvel2)>=0.0)
    wdir=1.0;

	 
	 d->M.p[count] = udir*ivel2/p->DXN[IM1] - (1.0-udir)*ivel1/p->DXN[IP]
					+ (vdir*jvel2/p->DYN[JM1] - (1.0-vdir)*jvel1/p->DYN[JP])*p->y_dir;
					+ wdir*kvel2/p->DZN[KM1] - (1.0-wdir)*kvel1/p->DZN[KP];
	 
	 d->M.s[count] = -udir*ivel1/p->DXN[IM1];
	 d->M.n[count] =  (1.0-udir)*ivel2/p->DXN[IP];
	 
	 d->M.e[count] = -vdir*jvel1/p->DYN[JM1]*p->y_dir;
	 d->M.w[count] =  (1.0-vdir)*jvel2/p->DYN[JP]*p->y_dir;
	 
	 d->M.b[count] = -wdir*kvel1/p->DZN[KM1];
	 d->M.t[count] =  (1.0-wdir)*kvel2/p->DZN[KP];
     
	 ++count;
    }
}
