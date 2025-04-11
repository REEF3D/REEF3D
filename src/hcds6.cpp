/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
Author: Elyas Larkermani
--------------------------------------------------------------------*/

#include"hcds6.h"
#include"lexer.h"
#include"fdm.h"
#include"flux.h"

hcds6::hcds6(flux* _pflux)
{
    pflux = _pflux;
}

hcds6::~hcds6()
{
}

void hcds6::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{
    if(ipol==1)
    ULOOP
    a->F(i,j,k)+=aij(p,a,b,1,uvel,vvel,wvel,p->DXP,p->DYN,p->DZN,p->DXN,p->DYP,p->DZP);

    if(p->j_dir==1)
    if(ipol==2)
    VLOOP
    a->G(i,j,k)+=aij(p,a,b,2,uvel,vvel,wvel,p->DXN,p->DYP,p->DZN,p->DXP,p->DYN,p->DZP);

    if(ipol==3)
    WLOOP
    a->H(i,j,k)+=aij(p,a,b,3,uvel,vvel,wvel,p->DXN,p->DYN,p->DZP,p->DXP,p->DYP,p->DZN);

    if(ipol==4)
    LOOP
    a->L(i,j,k)+=aij(p,a,b,4,uvel,vvel,wvel,p->DXN,p->DYN,p->DZN,p->DXP,p->DYP,p->DZP);
    
    if(ipol==5)
    LOOP
    a->L(i,j,k)+=aij(p,a,b,5,uvel,vvel,wvel,p->DXN,p->DYN,p->DZN,p->DXP,p->DYP,p->DZP);

}

double hcds6::aij(lexer* p,fdm* a,field& b,int ipol, field& uvel, field& vvel, field& wvel, double *DX,double *DY, double *DZ, double *DXX,double *DYY, double *DZZ)
{				
		dx=dy=dz=0.0;


        pflux->u_flux(a,ipol,uvel,ivel1,ivel2);
        pflux->v_flux(a,ipol,vvel,jvel1,jvel2);
        pflux->w_flux(a,ipol,wvel,kvel1,kvel2);	



		dx = (ivel2*(111.0*b(i,j,k) + 111.0*b(i+1,j,k) - 24.0*b(i+2,j,k) + 3.0*b(i+3,j,k) - 24.0*b(i-1,j,k) + 3.0*b(i-2,j,k))  
		 -  ivel1* (111.0*b(i,j,k) + 111.0*b(i-1,j,k) - 24.0*b(i-2,j,k) + 3.0*b(i-3,j,k) - 24.0*b(i+1,j,k) + 3.0*b(i+2,j,k)))/(180.0*DX[IP]);

		
		if(p->j_dir==1)
		dy = (jvel2*(111.0*b(i,j,k) + 111.0*b(i,j+1,k) - 24.0*b(i,j+2,k) + 3.0*b(i,j+3,k) - 24.0*b(i,j-1,k) + 3.0*b(i,j-2,k))  
		 -  jvel1* (111.0*b(i,j,k) + 111.0*b(i,j-1,k) - 24.0*b(i,j-2,k) + 3.0*b(i,j-3,k) - 24.0*b(i,j+1,k) + 3.0*b(i,j+2,k)))/(180.0*DY[JP]);


		
		dz = (kvel2*(111.0*b(i,j,k) + 111.0*b(i,j,k+1) - 24.0*b(i,j,k+2) + 3.0*b(i,j,k+3) - 24.0*b(i,j,k-1) + 3.0*b(i,j,k-2))  
		 -  kvel1* (111.0*b(i,j,k) + 111.0*b(i,j,k-1) - 24.0*b(i,j,k-2) + 3.0*b(i,j,k-3) - 24.0*b(i,j,k+1) + 3.0*b(i,j,k+2)))/(180.0*DZ[KP]);
		
	
    		
		L = -dx-dy-dz;

		return L;
}

