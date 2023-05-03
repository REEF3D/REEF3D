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

#include"flux_face_FOU.h"
#include"lexer.h"
#include"fdm.h"

flux_face_FOU::flux_face_FOU(lexer *pp)
{
    p=pp;

}

flux_face_FOU::~flux_face_FOU()
{
}

void flux_face_FOU::u_flux(fdm *a, int ipol, field& uvel, double &uflux1, double &uflux2)
{
	if(ipol==1)
	{
        pip=1;
        
        if(p->flag1[UIm1JK]>0)
        {
        if(0.5*(uvel(i,j,k)+uvel(i-1,j,k)) >= 0.0)
        uflux1 = uvel(i-1,j,k);
        
        if(0.5*(uvel(i,j,k)+uvel(i-1,j,k)) < 0.0)
        uflux1 = uvel(i,j,k);
        }
        
        if(p->flag1[UIm1JK]<0)
        uflux1= 0.5*(uvel(i,j,k)+uvel(i-1,j,k));
        
        
        if(p->flag1[UIp1JK]>0)
        {
        if(0.5*(uvel(i,j,k)+uvel(i+1,j,k)) >= 0.0)
        uflux2 = uvel(i,j,k);
        
        if(0.5*(uvel(i,j,k)+uvel(i+1,j,k)) < 0.0)
        uflux2 = uvel(i+1,j,k);
        }
        
        if(p->flag1[UIp1JK]<0)
        uflux2= 0.5*(uvel(i,j,k)+uvel(i+1,j,k));

        pip=0;
	}

	if(ipol==2)
	{
	pip=1;
	uflux1= 0.5*(uvel(i-1,j,k)+uvel(i-1,j+1,k));
	uflux2= 0.5*(uvel(i,j,k)+uvel(i,j+1,k));
	pip=0;
	}

	if(ipol==3)
	{
	pip=1;
	uflux1= 0.5*(uvel(i-1,j,k)+uvel(i-1,j,k+1));
	uflux2= 0.5*(uvel(i,j,k)+uvel(i,j,k+1));
	pip=0;
	}

	if(ipol==4)
	{
	uflux1= uvel(i-1,j,k);
	uflux2= uvel(i,j,k);
	}

}

void flux_face_FOU::v_flux(fdm *a, int ipol, field& vvel, double &vflux1, double &vflux2)
{
	if(ipol==1)
	{
	pip=2;
	vflux1= 0.5*(vvel(i,j-1,k)+vvel(i+1,j-1,k));
	vflux2= 0.5*(vvel(i,j,k)+vvel(i+1,j,k));
	pip=0;
	}

	if(ipol==2)
	{
        pip=2;
        
        if(p->flag2[VIJm1K]>0)
        {
        if(0.5*(vvel(i,j,k)+vvel(i,j-1,k)) >= 0.0)
        vflux1 = vvel(i,j-1,k);
        
        if(0.5*(vvel(i,j,k)+vvel(i,j-1,k)) < 0.0)
        vflux1 = vvel(i,j,k);
        }
        
        if(p->flag2[VIJm1K]<0)
        vflux1= 0.5*(vvel(i,j,k)+vvel(i,j-1,k));
            
            
            
        if(p->flag2[VIJp1K]>0)
        {
        if(0.5*(vvel(i,j,k)+vvel(i,j+1,k)) >= 0.0)
        vflux2 = vvel(i,j,k);
        
        if(0.5*(vvel(i,j,k)+vvel(i,j+1,k)) < 0.0)
        vflux2 = vvel(i,j+1,k);
        }
        
        if(p->flag2[VIJp1K]<0)
        vflux2= 0.5*(vvel(i,j,k)+vvel(i,j+1,k));
        
        pip=0;
	}

	if(ipol==3)
	{
	pip=2;
	vflux1= 0.5*(vvel(i,j-1,k)+vvel(i,j-1,k+1));
	vflux2= 0.5*(vvel(i,j,k)+vvel(i,j,k+1));
	pip=0;
	}

	if(ipol==4)
	{
    pip=2;
	vflux1= vvel(i,j-1,k);
	vflux2= vvel(i,j,k);
    pip=0;
	}
}

void flux_face_FOU::w_flux(fdm *a, int ipol, field& wvel, double &wflux1, double &wflux2)
{

	if(ipol==1)
	{
	pip=3;
	wflux1= 0.5*(wvel(i,j,k-1)+wvel(i+1,j,k-1));
	wflux2= 0.5*(wvel(i,j,k)+wvel(i+1,j,k));
	pip=0;
	}

	if(ipol==2)
	{
	pip=3;
	wflux1= 0.5*(wvel(i,j,k-1)+wvel(i,j+1,k-1));
	wflux2= 0.5*(wvel(i,j,k)+wvel(i,j+1,k));
	pip=0;
	}

	if(ipol==3)
	{
        pip=3;
        
        if(p->flag3[WIJKm1]>0)
        {
        if(0.5*(wvel(i,j,k)+wvel(i,j,k-1))>=0.0)
        wflux1 = wvel(i,j,k-1);
        
        if(0.5*(wvel(i,j,k)+wvel(i,j,k-1))<0.0)
        wflux1 = wvel(i,j,k);
        }
        
        if(p->flag3[WIJKm1]<0)
        wflux1= 0.5*(wvel(i,j,k)+wvel(i,j,k-1));
        
        if(p->flag3[WIJKp1]>0)
        {
        if(0.5*(wvel(i,j,k)+wvel(i,j,k+1))>=0.0)
        wflux2 = wvel(i,j,k);
        
        if(0.5*(wvel(i,j,k)+wvel(i,j,k))<0.0)
        wflux2 = wvel(i,j,k+1);
        }
        
        if(p->flag3[WIJKp1]<0)
        wflux2= 0.5*(wvel(i,j,k)+wvel(i,j,k+1));

        pip=0;
	}

	if(ipol==4)
	{
        pip=3;
        wflux1= wvel(i,j,k-1);
        wflux2= wvel(i,j,k);
        pip=0;
	}
}

void flux_face_FOU::omega_flux(lexer *p, fdm* a, int ipol, field &u, field &v, field& w, double &wflux1, double &wflux2)
{

	if(ipol==1)
	{
	pip=3;
	wflux1= 0.25*(w(i,j,k-1)+w(i+1,j,k-1)+w(i,j,k)+w(i+1,j,k));
	wflux2= 0.25*(w(i,j,k)+w(i+1,j,k)+w(i,j,k+1)+w(i+1,j,k+1));
	pip=0;
	}

	if(ipol==2)
	{
	pip=3;
	wflux1= 0.25*(w(i,j,k-1)+w(i,j+1,k-1)+w(i,j,k)+w(i,j+1,k));
	wflux2= 0.25*(w(i,j,k)+w(i,j+1,k)+w(i,j,k+1)+w(i,j+1,k+1));
	pip=0;
	}

	if(ipol==3)
	{
    pip=3;
	wflux1= w(i,j,k);
	wflux2= w(i,j,k+1);
	pip=0;
	}

	if(ipol==4)
	{
    pip=3;
	wflux1= 0.5*(w(i,j,k-1)+w(i,j,k));
	wflux2= 0.5*(w(i,j,k)+w(i,j,k+1));
    pip=0;
	}
}
