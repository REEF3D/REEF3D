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

#include"iquick.h"
#include"lexer.h"
#include"fdm.h"
#include"flux_face_CDS2.h"
#include"flux_face_CDS2_vrans.h"
#include"flux_face_FOU.h"
#include"flux_face_FOU_vrans.h"
#include"flux_face_QOU.h"

iquick::iquick (lexer *p)
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
}

iquick::~iquick()
{
}

void iquick::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
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

void iquick::aij(lexer* p,fdm* a,field& b,field &F,int ipol, field& uvel, field& vvel, field& wvel)
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
	
	F(i,j,k) -=  ((ivel1*ul)/(8.0*p->dx)*b(i-2,j,k)
				-(ivel2*(1.0-ur))/(8.0*p->dx)*b(i+2,j,k)
				+(jvel1*vl)/(8.0*p->dx)*b(i,j-2,k)
				-(jvel2*(1.0-vr))/(8.0*p->dx)*b(i,j+2,k)
				+(kvel1*wl)/(8.0*p->dx)*b(i,j,k-2)
				-(kvel2*(1.0-wr))/(8.0*p->dx)*b(i,j,k+2))/p->dx;

	 
	 a->M.p[count] = ((6.0*ivel2*ur + 3.0*ivel2*(1.0-ur) - 3.0*ivel1*ul - 6.0*ivel1*(1.0-ul))
					+ (6.0*jvel2*vr + 3.0*jvel2*(1.0-vr) - 3.0*jvel1*vl - 6.0*jvel1*(1.0-vl))
					+ (6.0*kvel2*wr + 3.0*kvel2*(1.0-wr) - 3.0*kvel1*wl - 6.0*kvel1*(1.0-wl)))/(8.0*p->dx);
	 
	 a->M.s[count] = (-ivel2*ur - 6.0*ivel1*ul + ivel1*(1.0-ul) )/(8.0*p->dx);
	 a->M.n[count] = (3.0*ivel2*ur + 6.0*ivel2*(1.0-ur) - 3.0*ivel1*(1.0-ul) )/(8.0*p->dx);
	 
	 a->M.e[count] = (-jvel2*vr - 6.0*jvel1*vl + jvel1*(1.0-vl) )/(8.0*p->dx);
	 a->M.w[count] = (3.0*jvel2*vr + 6.0*jvel2*(1.0-vr) - 3.0*jvel1*(1.0-vl) )/(8.0*p->dx);
	 
	 a->M.b[count] = (-kvel2*wr - 6.0*kvel1*wl + kvel1*(1.0-wl) )/(8.0*p->dx);
	 a->M.t[count] = (3.0*kvel2*wr + 6.0*kvel2*(1.0-wr) - 3.0*kvel1*(1.0-wl) )/(8.0*p->dx);
	 
	++count;
}



