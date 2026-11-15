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

#include"sflow_ifou.h"
#include"lexer.h"
#include"fdm2D.h"
#include"sflow_flux_face_CDS.h"

sflow_ifou::sflow_ifou(lexer* p)
{
    pflux = new sflow_flux_face_CDS(p);
}

sflow_ifou::~sflow_ifou()
{
}

void sflow_ifou::start(lexer* p, fdm2D* b, slice& f, int ipol, slice& uvel, slice& vvel)
{
    count=0;
    
    if(ipol==1)
    SLICELOOP1
    aij(p,b,f,1,uvel,vvel);

    if(ipol==2)
    SLICELOOP2
    aij(p,b,f,2,uvel,vvel);
    
    if(ipol==4)
    SLICELOOP4
    aij(p,b,f,4,uvel,vvel);
}

// Implicit first-order upwind in advective form, u df/dx + v df/dy, as the upwind flux of each face minus f times
// the divergence of the face velocities:
//   M.p = (u1+ - u2-)/dx + (v1+ - v2-)/dy,  M.s = -u1+/dx,  M.n = u2-/dx,  M.e = -v1+/dy,  M.w = v2-/dy
// (u1, v1 the west / south face, u2, v2 the east / north face). The rows sum to zero and the off-diagonals are
// not positive: a constant stays constant and the solution stays positive wherever the depth-averaged velocity
// has a divergence (a sloping bed, waves, wetting fronts), where the flux form d(uf)/dx adds the source -f div(u).
void sflow_ifou::aij(lexer* p,fdm2D* b,slice& f,int ipol, slice& uvel, slice& vvel)
{
    pflux->u_flux(ipol,uvel,ivel1,ivel2);
    pflux->v_flux(ipol,vvel,jvel1,jvel2);

    const double up1 = MAX(ivel1,0.0);
    const double um2 = MIN(ivel2,0.0);
    const double vp1 = MAX(jvel1,0.0);
    const double vm2 = MIN(jvel2,0.0);

    // control-volume widths
    const double dx = (ipol==1) ? p->DXP[IP] : p->DXN[IP];
    const double dy = (ipol==2) ? p->DYP[JP] : p->DYN[JP];

    b->M.p[count] = (up1 - um2)/dx + (vp1 - vm2)/dy*p->y_dir;

    b->M.s[count] = -up1/dx;
    b->M.n[count] =  um2/dx;

    b->M.e[count] = -vp1/dy*p->y_dir;
    b->M.w[count] =  vm2/dy*p->y_dir;

    ++count;
}
