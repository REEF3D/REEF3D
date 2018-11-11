/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"ihires.h"
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
#include"flux_face_CDS2_vrans.h"
#include"flux_face_FOU.h"
#include"flux_face_FOU_vrans.h"
#include"flux_face_QOU.h"

ihires::ihires (lexer *p, int limiter) 
{
    if(p->B269==0)
    {
        if(p->D11==1)
        pflux = new flux_face_FOU(p);
        
        if(p->D11==2)
        pflux = new flux_face_CDS2(p);
        
        if(p->D11==3)
        pflux = new flux_face_QOU(p);
    }
    
    if(p->B269>=1)
    {
        if(p->D11==1)
        pflux = new flux_face_FOU_vrans(p);
        
        if(p->D11==2)
        pflux = new flux_face_CDS2_vrans(p);
        
        if(p->D11==3)
        pflux = new flux_face_FOU_vrans(p);
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

ihires::~ihires()
{
}

void ihires::start(lexer* p, fdm* a,  field& b, int ipol, field& uvel, field& vvel, field& wvel)
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

void ihires::aij(lexer* p,fdm* a,field& b,field& F,int ipol, field& uvel, field& vvel, field& wvel)
{
    ul=ur=vl=vr=wl=wr=0.0;
		
	pflux->u_flux(a,ipol,uvel,ivel1,ivel2);
    pflux->v_flux(a,ipol,vvel,jvel1,jvel2);
    pflux->w_flux(a,ipol,wvel,kvel1,kvel2);

	if(ivel1>=0.0)
	ul=1.0;

	if(ivel2>=0.0)
	ur=1.0;

	if(jvel1>=0.0)
	vl=1.0;

	if(jvel2>=0.0)
	vr=1.0;

	if(kvel1>=0.0)
	wl=1.0;

	if(kvel2>=0.0)
	wr=1.0;
		
	F(i,j,k) -= (-(ivel2*(1.0-ur)*0.5*plim->iphi(b,1,0,2,1))*b(i+2,j,k)
				 -(jvel2*(1.0-vr)*0.5*plim->jphi(b,1,0,2,1))*b(i,j+2,k)
				 -(kvel2*(1.0-wr)*0.5*plim->kphi(b,1,0,2,1))*b(i,j,k+2))/p->dx;		

	 
	 a->M.s[count] = -ivel1*ul*(1.0 - 0.5*plim->iphi(b,-1,-2,0,-1))/p->dx;
	 a->M.n[count] =  (ivel2*((ur*0.5*plim->iphi(b,0,-1,1,0)) + (1.0-ur)*(1.0 + 0.5*plim->iphi(b,1,0,2,1))) 
		                  + ivel1*(1.0-ul)*0.5*plim->iphi(b,0,-1,1,0))/p->dx;
	 
	 a->M.e[count] = -jvel1*vl*(1.0 - 0.5*plim->jphi(b,-1,-2,0,-1))/p->dx;
	 a->M.w[count] =  (jvel2*((vr*0.5*plim->jphi(b,0,-1,1,0)) + (1.0-vr)*(1.0 + 0.5*plim->jphi(b,1,0,2,1))) 
		                  + jvel1*(1.0-vl)*0.5*plim->jphi(b,0,-1,1,0))/p->dx;
	 
	 a->M.b[count] = -kvel1*wl*(1.0 - 0.5*plim->kphi(b,-1,-2,0,-1))/p->dx;
	 a->M.t[count] =  (kvel2*((wr*0.5*plim->kphi(b,0,-1,1,0)) + (1.0-wr)*(1.0 + 0.5*plim->kphi(b,1,0,2,1))) 
		                  + kvel1*(1.0-wl)*0.5*plim->kphi(b,0,-1,1,0))/p->dx;
	 
	 ++count; 
}
