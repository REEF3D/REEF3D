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

#include"cds4.h"
#include"lexer.h"
#include"fdm.h"
#include"flux_face_CDS2.h"
#include"flux_face_CDS4.h"
#include"flux_face_CDS2_vrans.h"
#include"flux_face_FOU.h"
#include"flux_face_FOU_vrans.h"
#include"flux_face_QOU.h"

cds4::cds4 (lexer *p)
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

cds4::~cds4()
{
}

void cds4::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{
    if(ipol==1)
    ULOOP
    a->F(i,j,k)+=aij(p,a,b,1,uvel,vvel,wvel);

    if(ipol==2)
    VLOOP
    a->G(i,j,k)+=aij(p,a,b,2,uvel,vvel,wvel);

    if(ipol==3)
    WLOOP
    a->H(i,j,k)+=aij(p,a,b,3,uvel,vvel,wvel);

    if(ipol==4)
    LOOP
    a->L(i,j,k)+=aij(p,a,b,4,uvel,vvel,wvel);

    if(ipol==5)
    LOOP
    a->L(i,j,k)+=aij(p,a,b,5,uvel,vvel,wvel);

}

double cds4::aij(lexer* p,fdm* a,field& b,int ipol, field& uvel, field& vvel, field& wvel)
{		
		dx=dy=dz=0.0;
		
		pflux->u_flux(a,ipol,uvel,ivel1,ivel2);
        pflux->v_flux(a,ipol,vvel,jvel1,jvel2);
        pflux->w_flux(a,ipol,wvel,kvel1,kvel2);
        
		
		dx= (ivel2*(27.0*b(i,j,k) + 27.0*b(i+1,j,k) - 3.0*b(i+2,j,k) - 3.0*b(i-1,j,k))  
		 -  ivel1* (27.0*b(i,j,k) + 27.0*b(i-1,j,k) - 3.0*b(i-2,j,k) - 3.0*b(i+1,j,k)))/(48.0*p->DXM);

		
		
		dy= (jvel2*(27.0*b(i,j,k) + 27.0*b(i,j+1,k) - 3.0*b(i,j+2,k) - 3.0*b(i,j-1,k))  
		 -  jvel1* (27.0*b(i,j,k) + 27.0*b(i,j-1,k) - 3.0*b(i,j-2,k) - 3.0*b(i,j+1,k)))/(48.0*p->DXM);


		
		dz= (kvel2*(27.0*b(i,j,k) + 27.0*b(i,j,k+1) - 3.0*b(i,j,k+2) - 3.0*b(i,j,k-1))  
		 -  kvel1* (27.0*b(i,j,k) + 27.0*b(i,j,k-1) - 3.0*b(i,j,k-2) - 3.0*b(i,j,k+1)))/(48.0*p->DXM);
		
		L = -dx-dy-dz;

		return L;
}

