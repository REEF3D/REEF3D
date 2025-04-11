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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"quick.h"
#include"lexer.h"
#include"fdm.h"
#include"flux.h"

quick::quick(flux* _pflux)
{
    pflux = _pflux;
}

quick::~quick()
{

}

void quick::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{
    if(ipol==1)
    ULOOP
    a->F(i,j,k)+=aij(p,a,b,1,uvel,vvel,wvel);
    
    if(p->j_dir==1)
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

double quick::aij(lexer* p,fdm* a,field& b,int ipol, field& uvel, field& vvel, field& wvel)
{

    ul=ur=vl=vr=wl=wr=dx=dy=dz=0.0;
		
		pflux->u_flux(a,ipol,uvel,ivel1,ivel2);
        pflux->v_flux(a,ipol,vvel,jvel1,jvel2);
        pflux->w_flux(a,ipol,wvel,kvel1,kvel2);
		
		if(ivel1>=0.0)
		ul=1.0;

		if(ivel2>=0.0)
		ur=1.0;

		dx = (ivel2*(ur*(3.0*b(i+1,j,k) + 6.0*b(i,j,k) - b(i-1,j,k))
              +(1.0-ur)*(3.0*b(i,j,k) + 6.0*b(i+1,j,k) - b(i+2,j,k)))

           -  ivel1*(ul*(3.0*b(i,j,k) + 6.0*b(i-1,j,k) - b(i-2,j,k))
              +(1.0-ul)*(3.0*b(i-1,j,k) + 6.0*b(i,j,k) - b(i+1,j,k))))/(8.0*p->DXM);



		if(p->j_dir==1)
        {
		if(jvel1>=0.0)
		vl=1.0;

		if(jvel2>=0.0)
		vr=1.0;

		dy = (jvel2*(vr*(3.0*b(i,j+1,k) + 6.0*b(i,j,k) - b(i,j-1,k))
              +(1.0-vr)*(3.0*b(i,j,k) + 6.0*b(i,j+1,k) - b(i,j+2,k)))

           -  jvel1*(vl*(3.0*b(i,j,k) + 6.0*b(i,j-1,k) - b(i,j-2,k))
              +(1.0-vl)*(3.0*b(i,j-1,k) + 6.0*b(i,j,k) - b(i,j+1,k))))/(8.0*p->DXM);
        }



		if(kvel1>=0.0)
		wl=1.0;

		if(kvel2>=0.0)
		wr=1.0;

		dz = (kvel2*(wr*(3.0*b(i,j,k+1) + 6.0*b(i,j,k) - b(i,j,k-1))
              +(1.0-wr)*(3.0*b(i,j,k) + 6.0*b(i,j,k+1) - b(i,j,k+2)))

           -  kvel1*(wl*(3.0*b(i,j,k) + 6.0*b(i,j,k-1) - b(i,j,k-2))
              +(1.0-wl)*(3.0*b(i,j,k-1) + 6.0*b(i,j,k) - b(i,j,k+1))))/(8.0*p->DXM);
		

		L = -dx-dy-dz;

		return L;
}

