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

#include"ifou.h"
#include"lexer.h"
#include"fdm.h"
#include"flux_face_CDS2.h"
#include"flux_face_CDS2_vrans.h"
#include"flux_face_FOU.h"
#include"flux_face_FOU_vrans.h"
#include"flux_face_CDS2_2D.h"
#include"flux_face_CDS2_vrans_2D.h"
#include"flux_face_FOU_2D.h"
#include"flux_face_FOU_vrans_2D.h"

ifou::ifou (lexer *p)
{
    if(p->j_dir==0)
    {
    if(p->B200==0)
    {
        if(p->D11==1)
        pflux = new flux_face_FOU_2D(p);
        
        if(p->D11==2)
        pflux = new flux_face_CDS2_2D(p);
    }
    
    if(p->B200>=1 || p->S10==2)
    {
        if(p->D11==1)
        pflux = new flux_face_FOU_vrans_2D(p);
        
        if(p->D11==2)
        pflux = new flux_face_CDS2_vrans_2D(p);
    }
    }
    
    if(p->j_dir==1)
    {
    if(p->B200==0)
    {
        if(p->D11==1)
        pflux = new flux_face_FOU(p);
        
        if(p->D11==2)
        pflux = new flux_face_CDS2(p);
    }
    
    if(p->B200>=1 || p->S10==2)
    {
        if(p->D11==1)
        pflux = new flux_face_FOU_vrans(p);
        
        if(p->D11==2)
        pflux = new flux_face_CDS2_vrans(p);
    }
    }
}

ifou::~ifou()
{
}

void ifou::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{
    count=0;

    if(ipol==1)
    ULOOP
    aij(p,a,b,1,uvel,vvel,wvel,p->DXP,p->DYN,p->DZN);
    
    if(p->j_dir==1)
    if(ipol==2)
    VLOOP
    aij(p,a,b,2,uvel,vvel,wvel,p->DXN,p->DYP,p->DZN);

    if(ipol==3)
    WLOOP
    aij(p,a,b,3,uvel,vvel,wvel,p->DXN,p->DYN,p->DZP);

    if(ipol==4)
    LOOP
    aij(p,a,b,4,uvel,vvel,wvel,p->DXN,p->DYN,p->DZN);

    if(ipol==5)
    LOOP
    aij(p,a,b,5,uvel,vvel,wvel,p->DXN,p->DYN,p->DZN);
}

// first-order upwind, implicit, conservative: (F_i+1/2 - F_i-1/2)/dx_i with the upwind value chosen per face,
// F_i+1/2 = max(u,0) b_i + min(u,0) b_i+1; dx_i is the width of the control volume (DXN for cell values, DXP for
// u in x, ...). The coefficients form an M-matrix (diagonal >= 0, off-diagonals <= 0) on any grid.
// (Before: one upwind direction per cell from the mean of the two face velocities and, for positive flow, the
// width of the upstream cell DX[IM1], which is not conservative on stretched grids.)
void ifou::aij(lexer* p,fdm* a,field& b,int ipol, field& uvel, field& vvel, field& wvel, double *DX,double *DY, double *DZ)
{
    pflux->u_flux(a,ipol,uvel,ivel1,ivel2);
    pflux->v_flux(a,ipol,vvel,jvel1,jvel2);
    pflux->w_flux(a,ipol,wvel,kvel1,kvel2);
    
    const double dx = DX[IP];
    const double dy = DY[JP];
    const double dz = DZ[KP];

	 a->M.p[count] =  (MAX(ivel2,0.0) - MIN(ivel1,0.0))/dx
					+ (MAX(jvel2,0.0) - MIN(jvel1,0.0))/dy*p->y_dir
					+ (MAX(kvel2,0.0) - MIN(kvel1,0.0))/dz;
	 
	 a->M.s[count] = -MAX(ivel1,0.0)/dx;
	 a->M.n[count] =  MIN(ivel2,0.0)/dx;
	 
	 a->M.e[count] = -MAX(jvel1,0.0)/dy*p->y_dir;
	 a->M.w[count] =  MIN(jvel2,0.0)/dy*p->y_dir;
	 
	 a->M.b[count] = -MAX(kvel1,0.0)/dz;
	 a->M.t[count] =  MIN(kvel2,0.0)/dz;
     
	 ++count;
}
