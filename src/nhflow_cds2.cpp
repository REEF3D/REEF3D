/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"nhflow_cds2.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"nhflow_flux_face_cds2.h"

nhflow_cds2::nhflow_cds2 (lexer *p)
{
    pflux = new nhflow_flux_face_cds2(p);
}

nhflow_cds2::~nhflow_cds2()
{
}

void nhflow_cds2::start(lexer* p, fdm_nhf* d, double *F, int ipol, double *U, double *V, double *W)
{
    if(ipol==1)
        LOOP
        d->F[IJK]+=aij(p,d,F,1,U,V,W,p->DXN,p->DYN,p->DZN);

        if(ipol==2)
        LOOP
        d->G[IJK]+=aij(p,d,F,2,U,V,W,p->DXN,p->DYN,p->DZN);

        if(ipol==3)
        LOOP
        d->H[IJK]+=aij(p,d,F,3,U,V,W,p->DXN,p->DYN,p->DZN);

        if(ipol==4)
        LOOP
        d->L[IJK]+=aij(p,d,F,4,U,V,W,p->DXN,p->DYN,p->DZN);
}

double nhflow_cds2::aij(lexer* p,fdm_nhf* d, double *F, int ipol, double *U, double *V, double *W, double *DX,double *DY, double *DZ)
{
    udir=vdir=wdir=0.0;
    
        // convective flux
        pflux->u_flux(d,ipol,U,ivel1,ivel2);
        pflux->v_flux(d,ipol,V,jvel1,jvel2);
        pflux->w_flux(d,ipol,d->omegaF,kvel1,kvel2);
    
        dx = (ivel2*0.5*(F[IJK] + F[Ip1JK])  -  ivel1*0.5*(F[Im1JK]  +  F[IJK]))/DX[IP];
		
		
		dy = (jvel2*0.5*(F[IJK] + F[IJp1K])  -  jvel1*0.5*(F[IJm1K]  +  F[IJK]))/DY[JP];
		
	
		dz = (kvel2*0.5*(F[IJK] + F[IJKp1])  -  kvel1*0.5*(F[IJKm1] +  F[IJK]))/DZ[KP];
    
    
    L = -dx-dy-dz;
    
    return L;
    
}

