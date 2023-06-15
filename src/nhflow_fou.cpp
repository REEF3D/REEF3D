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

#include"nhflow_fou.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"nhflow_flux_face_cds2.h"

nhflow_fou::nhflow_fou (lexer *p)
{
    pflux = new nhflow_flux_face_cds2(p);
}

nhflow_fou::~nhflow_fou()
{
}

void nhflow_fou::start(lexer* p, fdm_nhf* d, double *F, int ipol, double *U, double *V, double *W, slice &eta)
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

double nhflow_fou::aij(lexer* p,fdm_nhf* d, double *F, int ipol, double *U, double *V, double *W, double *DX,double *DY, double *DZ)
{
    udir=vdir=wdir=0.0;
    
        // convective flux
        pflux->u_flux(d,ipol,U,ivel1,ivel2);
        pflux->v_flux(d,ipol,V,jvel1,jvel2);
        pflux->w_flux(d,ipol,d->omegaF,kvel1,kvel2);

    
    // x-dir
    if(0.5*(ivel1+ivel2)>=0.0)
    udir=1.0;
    
    dx =     udir*(ivel2*F[IJK] - ivel1*F[Im1JK])/DX[IM1] 
    
    +   (1.0-udir)*(ivel2*F[Ip1JK] - ivel1*F[IJK])/DX[IP]; 
    
    
    // y-dir
    if(0.5*(jvel1+jvel2)>=0.0)
    vdir=1.0;
    
    dy =     vdir*(jvel2*F[IJK] - jvel1*F[IJm1K])/DY[JM1] 
    
    +   (1.0-vdir)*(jvel2*F[IJp1K] - jvel1*F[IJK])/DY[JP]; 
    
    
    // z-dir
    if(0.5*(kvel1+kvel2)>=0.0)
    wdir=1.0;
    
    dz =     wdir*(kvel2*F[IJK] - kvel1*F[IJKm1])/DZ[KM1] 
    
    +   (1.0-wdir)*(kvel2*F[IJKp1] - kvel1*F[IJK])/DZ[KP]; 
    
    
    L = -dx-dy-dz;
    
    return L;
}

