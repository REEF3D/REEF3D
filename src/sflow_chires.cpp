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
--------------------------------------------------------------------*/

#include"sflow_chires.h"
#include"lexer.h"
#include"fdm2D.h"
#include"sflow_fluxlim_minmod.h"
#include"sflow_fluxlim_vanleer.h"
#include"sflow_fluxlim_smart.h"
#include"sflow_flux_face_FOU.h"
#include"sflow_flux_face_CDS.h"
#include"sflow_flux_HJ_CDS.h"
#include"sflow_flux_face_C_FOU.h"
#include"sflow_flux_face_C_CDS.h"
#include"sflow_flux_face_C_HJ.h"

#define HXIJ (fabs(b->hx(i,j))>1.0e-20?b->hx(i,j):1.0e20)
#define HYIJ (fabs(b->hy(i,j))>1.0e-20?b->hy(i,j):1.0e20)

#define HP (fabs(b->hp(i,j))>1.0e-20?b->hp(i,j):1.0e20)
#define HPI (fabs(b->hp(i+1,j))>1.0e-20?b->hp(i+1,j):1.0e20)
#define HPJ (fabs(b->hp(i,j+1))>1.0e-20?b->hp(i,j+1):1.0e20)

#define HPX (0.5*(HP + HPI))
#define HPY (0.5*(HP + HPJ))
 

sflow_chires::sflow_chires (lexer *p, fdm2D *b, int limiter)
{
	if(limiter==6)
	plim = new sflow_fluxlim_minmod(p);
    
    if(limiter==7)
	plim = new sflow_fluxlim_vanleer(p);
    
    if(limiter==8)
	plim = new sflow_fluxlim_smart(p);
    
    if(p->A216==1)
    pflux = new sflow_flux_face_C_FOU(p,b);
        
    if(p->A216==2)
    pflux = new sflow_flux_face_C_CDS(p,b);
        
    if(p->A216==4)
    pflux = new sflow_flux_face_C_HJ(p,b);
    
	/*
	if(limiter==17)
	plim = new vanleer(p);
	
	if(limiter==18)
	plim = new smart(p);
	*/
}

sflow_chires::~sflow_chires()
{

}

void sflow_chires::start(lexer* p, fdm2D* b, slice& f, int ipol, slice& uvel, slice& vvel)
{ 	
    if(ipol==1)
    SLICELOOP1
    b->F(i,j)+=aij(p,b,f,1,uvel,vvel);

    if(ipol==2)
    SLICELOOP2
    b->G(i,j)+=aij(p,b,f,2,uvel,vvel);
    
    if(ipol==4)
    SLICELOOP4
    b->L(i,j)+=aij(p,b,f,4,uvel,vvel);

    if(ipol==5)
    SLICELOOP4
    b->L(i,j)+=aij(p,b,f,5,uvel,vvel);

}

double sflow_chires::aij(lexer* p,fdm2D* b,slice& f,int ipol, slice& uvel, slice& vvel)
{

    ul=ur=vl=vr=wl=wr=dx=dy=0.0;
		
        pflux->u_flux(ipol,uvel,ivel1,ivel2);
        pflux->v_flux(ipol,vvel,jvel1,jvel2);

		
		if(ivel1>=0.0)
		ul=1.0;

		if(ivel2>=0.0)
		ur=1.0;

		dx = (ivel2*(ur*(f(i,j) + 0.5*plim->iphi(f,0,-1,1,0)*(f(i+1,j)-f(i,j)))
             +(1.0-ur)*(f(i+1,j) - 0.5*plim->iphi(f,1,0,2,1)*(f(i+2,j)-f(i+1,j))))

          -  ivel1*(ul*(f(i-1,j) + 0.5*plim->iphi(f,-1,-2,0,-1)*(f(i,j)-f(i-1,j)))
             +(1.0-ul)*(f(i,j) - 0.5*plim->iphi(f,0,-1,1,0)*(f(i+1,j)-f(i,j)))))/(p->DXM);
             
        dx -= f(i,j)*(ivel2-ivel1)/p->DXM;
        
            if(ipol==1)
            dx/=HXIJ;
            
            if(ipol==2)
            dx/=HYIJ;
            
            if(ipol==4)
            dx/=HP;



		if(jvel1>=0.0)
		vl=1.0;

		if(jvel2>=0.0)
		vr=1.0;

		dy = (jvel2*(vr*(f(i,j) + 0.5*plim->jphi(f,0,-1,1,0)*(f(i,j+1)-f(i,j)))
             +(1.0-vr)*(f(i,j+1) - 0.5*plim->jphi(f,1,0,2,1)*(f(i,j+2)-f(i,j+1))))

          -  jvel1*(vl*(f(i,j-1) + 0.5*plim->jphi(f,-1,-2,0,-1)*(f(i,j)-f(i,j-1)))
             +(1.0-vl)*(f(i,j) - 0.5*plim->jphi(f,0,-1,1,0)*(f(i,j)-f(i+1,j)))))/(p->DXM);
        
        dy -= f(i,j)*(jvel2-jvel1)/p->DXM;
        
            if(ipol==1)
            dy/=HXIJ;
            
            if(ipol==2)
            dy/=HYIJ;
        
            if(ipol==4)
            dy/=HP;

		
		L = -dx-dy;

		return L;
}
