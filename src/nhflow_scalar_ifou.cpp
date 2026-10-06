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
    
    // conservative form: (F_i+1/2 - F_i-1/2)/dx_i with the upwind value per face (as ifou); was one direction
    // per cell and the width of the upstream cell for positive flow (not conservative on stretched grids)
    count=0;
    LOOP
    {
    padvec->uadvec(ipol,U,ivel1,ivel2);
    padvec->vadvec(ipol,V,jvel1,jvel2);
    padvec->wadvec(ipol,W,kvel1,kvel2);
    
    const double dxc = p->DXN[IP];
    const double dyc = p->DYN[JP];
    const double dzc = p->DZN[KP]*(p->WL[IJ]>1.0e-20?p->WL[IJ]:1.0e20);

	 d->M.p[count] =    (MAX(ivel2,0.0) - MIN(ivel1,0.0))/dxc
					+ (MAX(jvel2,0.0) - MIN(jvel1,0.0))/dyc*p->y_dir
					+  (MAX(kvel2,0.0) - MIN(kvel1,0.0))/dzc;
	 
	 d->M.s[count] = -MAX(ivel1,0.0)/dxc;
	 d->M.n[count] =  MIN(ivel2,0.0)/dxc;
	 
	 d->M.e[count] = -MAX(jvel1,0.0)/dyc*p->y_dir;
	 d->M.w[count] =  MIN(jvel2,0.0)/dyc*p->y_dir;
	 
	 d->M.b[count] = -MAX(kvel1,0.0)/dzc;
	 d->M.t[count] =  MIN(kvel2,0.0)/dzc;
     
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
