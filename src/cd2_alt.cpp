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

#include"cds2_alt.h"
#include"lexer.h"
#include"fdm.h"
#include"flux_HJ_CDS2.h"
#include"flux_HJ_CDS4.h"
#include"flux_HJ_CDS2_vrans.h"

cds2_alt::cds2_alt (lexer *p)
{
    if(p->B269==0 && p->D11!=4)
    pflux = new flux_HJ_CDS2(p);
    
    if(p->B269==0 && p->D11==4)
    pflux = new flux_HJ_CDS2(p);
    
    if(p->B269>=1 || p->S10==2)
    pflux = new flux_HJ_CDS2_vrans(p);
}

cds2_alt::~cds2_alt()
{
}

void cds2_alt::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{
    if(ipol==1)
    ULOOP
    a->F(i,j,k)+=aij(p,a,b,1,uvel,vvel,wvel,p->DXP,p->DYN,p->DZN);

    if(ipol==2)
    VLOOP
    a->G(i,j,k)+=aij(p,a,b,2,uvel,vvel,wvel,p->DXN,p->DYP,p->DZN);

    if(ipol==3)
    WLOOP
    a->H(i,j,k)+=aij(p,a,b,3,uvel,vvel,wvel,p->DXN,p->DYN,p->DZP);

    if(ipol==4)
    LOOP
    a->L(i,j,k)+=aij(p,a,b,4,uvel,vvel,wvel,p->DXN,p->DYN,p->DZN);
    
    if(ipol==5)
    LOOP
    a->L(i,j,k)+=aij(p,a,b,5,uvel,vvel,wvel,p->DXN,p->DYN,p->DZN);

}

double cds2_alt::aij(lexer* p,fdm* a,field& b,int ipol, field& uvel, field& vvel, field& wvel, double *DX,double *DY, double *DZ)
{		
		dx=dy=dz=0.0;	

		pflux->u_flux(a,ipol,uvel,iadvec,ivel2);
        pflux->v_flux(a,ipol,vvel,jadvec,jvel2);
        pflux->w_flux(a,ipol,wvel,kadvec,kvel2);
		
        
		dx = iadvec*(0.5*(b(i,j,k) + b(i+1,j,k))  -  0.5*(b(i-1,j,k) +  b(i,j,k)))/DX[IP];
		
		dy = jadvec*(0.5*(b(i,j,k) + b(i,j+1,k))  -  0.5*(b(i,j-1,k) +  b(i,j,k)))/DY[JP];
	        
		dz = kadvec*(0.5*(b(i,j,k) + b(i,j,k+1))  -  0.5*(b(i,j,k-1) +  b(i,j,k)))/DZ[KP];

		L = -dx-dy-dz;

		return L;
}

