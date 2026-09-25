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

#include"6DOF_obj.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void sixdof_obj::hydrodynamic_viscous_forces_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL, 
                                            double &Fv_x, double &Fv_y, double &Fv_z, double A_triang,
                                            double xp, double yp, double zp, double nx, double ny, double nz)
{
    // Simple skin-friction estimate on one surface element of area A_triang:
    // F = 0.5*rho*Cf*|u_t|*u_t*A, with u_t the fluid velocity relative to the local rigid-body
    // velocity, projected onto the element plane.
    // The outputs are always defined on return (X38!=1 -> no viscous load).
    Fv_x = Fv_y = Fv_z = 0.0;
    
    if(p->X38!=1)
    return;
    
    const double Cf=0.0001;
    
    // cell containing the surface point (the member i,j,k are stale here)
    const int ii = MAX(0, MIN(p->knox-1, p->posc_i(xp)));
    const int jj = (p->j_dir==1) ? MAX(0, MIN(p->knoy-1, p->posc_j(yp))) : 0;
    const int kk = MAX(0, MIN(p->knoz-1, p->posc_sig(ii,jj,zp)));
    
    // sample point X43 cells off the wall along the outward normal
    const double xs = xp + p->X43*nx*p->DXP[ii+marge];
    const double ys = yp + p->X43*ny*p->DYP[jj+marge];
    const double zs = zp + p->X43*nz*p->DZP[kk+marge]*WL(ii,jj);
    
    const double uval = p->ccipol4V(d->U, WL, d->bed, xs, ys, zs); 
    const double vval = p->ccipol4V(d->V, WL, d->bed, xs, ys, zs); 
    const double wval = p->ccipol4V(d->W, WL, d->bed, xs, ys, zs); 
    
    // rigid-body velocity at the sample point
    const double rx = xs - c_(0);
    const double ry = ys - c_(1);
    const double rz = zs - c_(2);
    
    double ur = uval - (u_fb(0) + u_fb(4)*rz - u_fb(5)*ry);
    double vr = vval - (u_fb(1) + u_fb(5)*rx - u_fb(3)*rz);
    double wr = wval - (u_fb(2) + u_fb(3)*ry - u_fb(4)*rx);
    
    // tangential part
    const double un = ur*nx + vr*ny + wr*nz;
    ur -= un*nx;
    vr -= un*ny;
    wr -= un*nz;
    
    const double ut = sqrt(ur*ur + vr*vr + wr*wr);
    
	Fv_x = 0.5*p->W1*Cf*ut*ur*A_triang;
    Fv_y = 0.5*p->W1*Cf*ut*vr*A_triang;
    Fv_z = 0.5*p->W1*Cf*ut*wr*A_triang;
    
    if(p->j_dir==0)
    Fv_y = 0.0;
}
