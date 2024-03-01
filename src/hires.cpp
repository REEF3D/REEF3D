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

#include"hires.h"
#include"lexer.h"
#include"fdm.h"
#include"minmod.h"
#include"vanleer.h"
#include"umist.h"
#include"vanalbada.h"
#include"superbee.h"
#include"smart.h"
#include"limo3.h"
#include"tvdvof.h"
#include"flux_face_CDS2.h"
#include"flux_face_CDS4.h"
#include"flux_face_CDS2_vrans.h"
#include"flux_face_FOU.h"
#include"flux_face_FOU_vrans.h"
#include"flux_face_QOU.h"

hires::hires (lexer *p, int limiter) 
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
    
    
	if(limiter==10)
	plim = new minmod(p);
	
	if(limiter==11)
	plim = new vanleer(p);
	
	if(limiter==12)
	plim = new umist(p);
	
	if(limiter==13)
	plim = new vanalbada(p);
	
	if(limiter==14)
	plim = new superbee(p);
	
	if(limiter==15)
	plim = new smart(p);
	
	if(limiter==16)
	plim = new limo3(p);
	
	if(limiter==42)
	plim = new tvdvof(p);
}

hires::~hires()
{

}

void hires::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
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

double hires::aij(lexer* p,fdm* a,field& b,int ipol, field& uvel, field& vvel, field& wvel, double *DX,double *DY, double *DZ)
{

    udir=vdir=wdir=0.0;

    pflux->u_flux(a,ipol,uvel,ivel1,ivel2);
    pflux->v_flux(a,ipol,vvel,jvel1,jvel2);
    pflux->w_flux(a,ipol,wvel,kvel1,kvel2);
		
		// x-dir
        if(0.5*(ivel1+ivel2)>=0.0)
        udir=1.0;

		dx = udir*(ivel2*(b(i,j,k) + 0.5*plim->iphi(b,0,-1,1,0)*(b(i+1,j,k)-b(i,j,k)))
        
                - ivel1*(b(i-1,j,k) + 0.5*plim->iphi(b,-1,-2,0,-1)*(b(i,j,k)-b(i-1,j,k))))/DX[IM1] 
            
            
             + (1.0-udir)*(ivel2*(b(i+1,j,k) - 0.5*plim->iphi(b,1,0,2,1)*(b(i+2,j,k)-b(i+1,j,k)))
          
             -            ivel1*(b(i,j,k) - 0.5*plim->iphi(b,0,-1,1,0)*(b(i+1,j,k)-b(i,j,k))))/DX[IP]; 
             

		// y-dir
        if(0.5*(jvel1+jvel2)>=0.0)
        vdir=1.0;

		dy = vdir*(jvel2*(b(i,j,k) + 0.5*plim->jphi(b,0,-1,1,0)*(b(i,j+1,k)-b(i,j,k)))
        
                - jvel1*(b(i,j-1,k) + 0.5*plim->jphi(b,-1,-2,0,-1)*(b(i,j,k)-b(i,j-1,k))))/DY[JM1] 
        
        
             + (1.0-vdir)*(jvel2*(b(i,j+1,k) - 0.5*plim->jphi(b,1,0,2,1)*(b(i,j+2,k)-b(i,j+1,k)))

             -             jvel1*(b(i,j,k) - 0.5*plim->jphi(b,0,-1,1,0)*(b(i,j,k)-b(i+1,j,k))))/DY[JP]; 


		// z-dir
        if(0.5*(kvel1+kvel2)>=0.0)
        wdir=1.0;

		dz = wdir*(kvel2*(b(i,j,k) + 0.5*plim->kphi(b,0,-1,1,0)*(b(i,j,k+1)-b(i,j,k)))
        
                -  kvel1*(b(i,j,k-1) + 0.5*plim->kphi(b,-1,-2,0,-1)*(b(i,j,k)-b(i,j,k-1))))/DZ[KM1] 
        
        
            + (1.0-wdir)*(kvel2*(b(i,j,k+1) - 0.5*plim->kphi(b,1,0,2,1)*(b(i,j,k+2)-b(i,j,k+1)))
          
             -      kvel1*(b(i,j,k) - 0.5*plim->kphi(b,0,-1,1,0)*(b(i,j,k+1)-b(i,j,k))))/DZ[KP];
		

		L = -dx-dy-dz;

		return L;
}
