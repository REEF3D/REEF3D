/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
#include"flux_face_CDS4.h"
#include"flux_face_CDS2_vrans.h"
#include"flux_face_FOU.h"
#include"flux_face_FOU_vrans.h"
#include"flux_face_QOU.h"

ifou::ifou (lexer *p)
{
    if(p->B269==0)
    {
        if(p->D11==1)
        pflux = new flux_face_FOU(p);
        
        if(p->D11==2)
        pflux = new flux_face_CDS2(p);
        
        if(p->D11==3)
        pflux = new flux_face_QOU(p);
        
        if(p->D11==4)
        pflux = new flux_face_CDS4(p);
    }
    
    if(p->B269>=1 || p->S10==2)
    {
        if(p->D11==1)
        pflux = new flux_face_FOU_vrans(p);
        
        if(p->D11==2)
        pflux = new flux_face_CDS2_vrans(p);
        
        if(p->D11==3)
        pflux = new flux_face_FOU_vrans(p);
        
        if(p->D11==4)
        pflux = new flux_face_CDS2(p);
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

void ifou::aij(lexer* p,fdm* a,field& b,int ipol, field& uvel, field& vvel, field& wvel, double *DX,double *DY, double *DZ)
{
	udir=vdir=wdir=0.0;
    
    pflux->u_flux(a,ipol,uvel,ivel1,ivel2);
    pflux->v_flux(a,ipol,vvel,jvel1,jvel2);
    pflux->w_flux(a,ipol,wvel,kvel1,kvel2);

	if(0.5*(ivel1+ivel2)>=0.0)
    udir=1.0;
    
    if(0.5*(jvel1+jvel2)>=0.0)
    vdir=1.0;
    
    if(0.5*(kvel1+kvel2)>=0.0)
    wdir=1.0;

	 
	 a->M.p[count] = udir*ivel2/DX[IM1] - (1.0-udir)*ivel1/DX[IP]
					+ (vdir*jvel2/DY[JM1] - (1.0-vdir)*jvel1/DY[JP])*p->y_dir;
					+ wdir*kvel2/DZ[KM1] - (1.0-wdir)*kvel1/DZ[KP];
	 
	 a->M.s[count] = -udir*ivel1/DX[IM1];
	 a->M.n[count] =  (1.0-udir)*ivel2/DX[IP];
	 
	 a->M.e[count] = -vdir*jvel1/DY[JM1]*p->y_dir;
	 a->M.w[count] =  (1.0-vdir)*jvel2/DY[JP]*p->y_dir;
	 
	 a->M.b[count] = -wdir*kvel1/DZ[KM1];
	 a->M.t[count] =  (1.0-wdir)*kvel2/DZ[KP];
     
	 ++count;
}




