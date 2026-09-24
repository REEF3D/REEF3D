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

#include"nhflow_scalar_ifou.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"nhflow_scalar_advec_CDS2.h"

nhflow_scalar_ifou::nhflow_scalar_ifou(lexer *p, int form) : advective(form)
{

    padvec = new nhflow_scalar_advec_CDS2(p);

}

nhflow_scalar_ifou::~nhflow_scalar_ifou()
{
}

void nhflow_scalar_ifou::start(lexer* p, fdm_nhf *d, double *F, int ipol, double *U, double *V, double *W)
{
    if(advective==1)
    {
    start_advective(p,d,F,ipol,U,V,W);
    return;
    }
    
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

	 
	 d->M.p[count] =    udir*ivel2/p->DXN[IM1] - (1.0-udir)*ivel1/p->DXN[IP]
					+ (vdir*jvel2/p->DYN[JM1] - (1.0-vdir)*jvel1/p->DYN[JP])*p->y_dir
					+  wdir*kvel2/(p->DZN[KM1]*p->WL[IJ]) - (1.0-wdir)*kvel1/(p->DZN[KP]*p->WL[IJ]);
	 
	 d->M.s[count] = -udir*ivel1/p->DXN[IM1];
	 d->M.n[count] =  (1.0-udir)*ivel2/p->DXN[IP];
	 
	 d->M.e[count] = -vdir*jvel1/p->DYN[JM1]*p->y_dir;
	 d->M.w[count] =  (1.0-vdir)*jvel2/p->DYN[JP]*p->y_dir;
	 
	 d->M.b[count] = -wdir*kvel1/(p->DZN[KM1]*p->WL[IJ]);
	 d->M.t[count] =  (1.0-wdir)*kvel2/(p->DZN[KP]*p->WL[IJ]);
     
	 ++count;
    }
}

void nhflow_scalar_ifou::start_advective(lexer* p, fdm_nhf *d, double *F, int ipol, double *U, double *V, double *W)
{
    // Implicit first-order upwind for the non-conservative (in D) scalar equation in sigma coordinates:
    //   dF/dt + u dF/dx|s + v dF/dy|s + (omega/D) dF/ds = ...
    // Written as per-face upwind flux divergence minus F*div(face velocities):
    //   M.p = (u1+ - u2-)/dx + ...,  M.s = -u1+/dx,  M.n = u2-/dx
    // -> row sum zero (constants preserved), M-matrix (positivity), no spurious F*dD/dt source.
    double up1,um2,vp1,vm2,wp1,wm2;
    double dxc,dyc,dzc;
    
    count=0;
    LOOP
    {
    padvec->uadvec(ipol,U,ivel1,ivel2);
    padvec->vadvec(ipol,V,jvel1,jvel2);
    padvec->wadvec(ipol,W,kvel1,kvel2);
    
    up1 = MAX(ivel1,0.0);   // inflow through i-1/2
    um2 = MIN(ivel2,0.0);   // inflow through i+1/2
    vp1 = MAX(jvel1,0.0);
    vm2 = MIN(jvel2,0.0);
    wp1 = MAX(kvel1,0.0);
    wm2 = MIN(kvel2,0.0);
    
    dxc = p->DXN[IP];
    dyc = p->DYN[JP];
    dzc = p->DZN[KP]*(p->WL[IJ]>1.0e-20?p->WL[IJ]:1.0e20);
    
	 d->M.p[count] =    (up1 - um2)/dxc
					+  ((vp1 - vm2)/dyc)*p->y_dir
					+   (wp1 - wm2)/dzc;
	 
	 d->M.s[count] = -up1/dxc;
	 d->M.n[count] =  um2/dxc;
	 
	 d->M.e[count] = -vp1/dyc*p->y_dir;
	 d->M.w[count] =  vm2/dyc*p->y_dir;
	 
	 d->M.b[count] = -wp1/dzc;
	 d->M.t[count] =  wm2/dzc;
     
	 ++count;
    }
}
