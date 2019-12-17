/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"icds4.h"
#include"lexer.h"
#include"fdm.h"
#include"flux_face_CDS2.h"
#include"flux_face_CDS4.h"
#include"flux_face_CDS2_vrans.h"
#include"flux_face_FOU.h"
#include"flux_face_FOU_vrans.h"
#include"flux_face_QOU.h"

icds4::icds4 (lexer *p)
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

icds4::~icds4()
{
}

void icds4::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{
    count=0;

    if(ipol==1)
    ULOOP
    aij(p,a,b,a->F,1,uvel,vvel,wvel);

    if(ipol==2)
    VLOOP
    aij(p,a,b,a->G,2,uvel,vvel,wvel);

    if(ipol==3)
    WLOOP
    aij(p,a,b,a->H,3,uvel,vvel,wvel);

    if(ipol==4)
    LOOP
    aij(p,a,b,a->L,4,uvel,vvel,wvel);

    if(ipol==5)
    LOOP
    aij(p,a,b,a->L,5,uvel,vvel,wvel);
}

void icds4::aij(lexer* p,fdm* a,field& b,field& F,int ipol, field& uvel, field& vvel, field& wvel)
{		
 	
    pflux->u_flux(a,ipol,uvel,ivel1,ivel2);
    pflux->v_flux(a,ipol,vvel,jvel1,jvel2);
    pflux->w_flux(a,ipol,wvel,kvel1,kvel2);
	
	F(i,j,k) -=  (3.0*ivel1)/(p->dx*48.0)*b(i-2,j,k)
				-(3.0*ivel2)/(p->dx*48.0)*b(i+2,j,k)
				+(3.0*jvel1)/(p->dx*48.0)*b(i,j-2,k)
				-(3.0*jvel2)/(p->dx*48.0)*b(i,j+2,k)
				+(3.0*kvel1)/(p->dx*48.0)*b(i,j,k-2)
				-(3.0*kvel2)/(p->dx*48.0)*b(i,j,k+2);
	 
	 
	 a->M.p[count] =    (27.0*ivel2 - 27.0*ivel1
						+ 27.0*jvel2 - 27.0*jvel1
						+ 27.0*kvel2 - 27.0*kvel1)/(p->dx*48.0);
	 
	 a->M.s[count] = (-3.0*ivel2 - 27.0*ivel1)/(p->dx*48.0);
	 a->M.n[count] = (27.0*ivel2 + 3.0*ivel1)/(p->dx*48.0);
	 
	 a->M.e[count] = (-3.0*jvel2 - 27.0*jvel1)/(p->dx*48.0);
	 a->M.w[count] = (27.0*jvel2 + 3.0*jvel1)/(p->dx*48.0);
	 
	 a->M.b[count] = (-3.0*kvel2 - 27.0*kvel1)/(p->dx*48.0);
	 a->M.t[count] = (27.0*kvel2 + 3.0*kvel1)/(p->dx*48.0);
	
	++count;
}





