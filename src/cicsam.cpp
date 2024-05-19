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

#include"cicsam.h"
#include"lexer.h"
#include"fdm.h"
#include"flux_face_CDS2.h"
#include"flux_face_CDS4.h"
#include"flux_face_CDS2_vrans.h"
#include"flux_face_FOU.h"
#include"flux_face_FOU_vrans.h"
#include"flux_face_QOU.h"

cicsam::cicsam (lexer *p) 
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

cicsam::~cicsam()
{
}

void cicsam::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
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

double cicsam::aij(lexer* p,fdm* a,field& b,int ipol, field& uvel, field& vvel, field& wvel)
{
    ul=ur=vl=vr=wl=wr=dx=dy=dz=0.0;
		

		pflux->u_flux(a,ipol,uvel,ivel1,ivel2);
        pflux->v_flux(a,ipol,vvel,jvel1,jvel2);
        pflux->w_flux(a,ipol,wvel,kvel1,kvel2);
        
        
        fx1 = cface(p,a,b,1,-1,ivel1);
        fx2 = cface(p,a,b,1,0,ivel2);

		dx= (ivel2*fx2 - ivel1*fx1)/(p->DXM);
		
		

		fy1 = cface(p,a,b,2,-1,jvel1);
	    fy2 = cface(p,a,b,2,0,jvel2);

		dy= (jvel2*fy2 - jvel1*fy1)/(p->DXM);
		
		
		
		fz1 = cface(p,a,b,3,-1,kvel1);
	    fz2 = cface(p,a,b,3,0,kvel2);

		dz= (kvel2*fz2 - kvel1*fz1)/(p->DXM);

		
		L = -dx-dy-dz;

		return L;
}

double cicsam::cface(lexer *p,fdm *a,field& b,int dir, int pos, double uwind)
{
	double cj,cj_,cj_s,cj_ss;
	double cc,cc_,cu,cd;
    double umax,Co,theta,costheta,gamma,beta;
	double gradx,grady,gradz;
	
	
	if(dir==1)
	{
		if(uwind>=0.0)
		{
		cc = b(i+pos,j,k);
		
		cu = b(i-1+pos,j,k);
		
		cd = b(i+1+pos,j,k);
        
        umax = 0.5*(a->u(i+pos,j,k)+a->u(i-1+pos,j,k));
		}
        
        if(uwind<0.0)
		{
		cc = b(i+1+pos,j,k);
		
		cu = b(i+2+pos,j,k);
		
		cd = b(i-0+pos,j,k);
        
        umax = 0.5*(a->u(i+1-pos,j,k)+a->u(i-0-pos,j,k));
		}
	}
	
	if(dir==2)
	{
		if(uwind>=0.0)
		{
		cc = b(i,j+pos,k);
		
		cu = b(i,j-1+pos,k);
		
		cd = b(i,j+1+pos,k);
        
        umax = 0.5*(a->v(i,j+pos,k)+a->v(i,j-1+pos,k));
		}
        
        if(uwind<0.0)
		{
		cc = b(i,j+1+pos,k);
	
		cu = b(i,j+2+pos,k);
		
		cd = b(i,j-0+pos,k);
        
        umax = 0.5*(a->v(i,j+1-pos,k)+a->v(i,j-0-pos,k));
		}
	}
	
	if(dir==3)
	{
		if(uwind>=0.0)
		{
		cc = b(i,j,k+pos);
		
		cu = b(i,j,k-1+pos);
		
		cd = b(i,j,k+1+pos);
        
        umax = 0.5*(a->w(i,j,k+pos)+a->w(i,j,k-1+pos));
		}
        
        if(uwind<0.0)
		{
		cc = b(i,j,k+1+pos);
		
		cu = b(i,j,k+2+pos);
		
		cd = b(i,j,k-0+pos);
        
        umax = 0.5*(a->w(i,j,k+1-pos)+a->w(i,j,k-0-pos));
		}
	}
	
	
	cc_ = (cc-cu)/(fabs(cd-cu)>1.0e-20?(cd-cu):1.0e20); 
    
    
    Co = fabs(umax*p->dt/p->DXM);
    
    
    if(cc_>=0.0 && cc_<1.0)
    cj_ = MIN(cc_/(fabs(Co)>1.0e-20?Co:1.0e20),1.0);
    
    if(cc_<0.0 || cc_>=1.0)
    cj_ = cc_;
    
    
    if(cc_>=0.0 && cc_<1.0)
    cj_s = MIN( (8.0*Co*cc_ + (1.0-Co)*(6.0*cc_+3.0))/8.0 , cj_);
    
    if(cc_<0.0 || cc_>=1.0)
    cj_s = cc_;
	
	
	
	if(uwind>=0.0)
	{
	gradx = fabs((b(i+1+pos,j,k)-b(i-1+pos,j,k))/(2.0*p->DXM));
	grady = fabs((b(i,j+1+pos,k)-b(i,j-1+pos,k))/(2.0*p->DXM));
	gradz = fabs((b(i,j,k+1+pos)-b(i,j,k-1+pos))/(2.0*p->DXM));
	}
	
	if(uwind<0.0)
	{
	gradx = fabs((b(i+1+1+pos,j,k)-b(i-1+1+pos,j,k))/(2.0*p->DXM));
	grady = fabs((b(i,j+1+1+pos,k)-b(i,j-1+1+pos,k))/(2.0*p->DXM));
	gradz = fabs((b(i,j,k+1+1+pos)-b(i,j,k-1+1+pos))/(2.0*p->DXM));
	}
	
	vl = sqrt(gradx*gradx + grady*grady + gradz*gradz);
    

	if(dir==1)
	theta = acos(gradx/(vl>1.0e-20?vl:1.0e20));

	if(dir==2)
	theta = acos(grady/(vl>1.0e-20?vl:1.0e20));
	
	if(dir==3)
	theta = acos(gradz/(vl>1.0e-20?vl:1.0e20));
	
	
    
    gamma = MIN((cos(2.0*theta)+1.0)/2.0,1.0);
    
    
	
    cj_ss = cj_*gamma + cj_s*(1.0-gamma);
    
    
    beta = (cj_ss-cc_)/(fabs(1.0-cc_)>1.0e-20?(1.0-cc_):1.0e20);
	
    
    cj = (1.0-beta)*cc + beta*cd;


    return cj;
}
