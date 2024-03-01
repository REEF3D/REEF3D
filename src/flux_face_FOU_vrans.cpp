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
#include"flux_face_FOU_vrans.h"
#include"lexer.h"
#include"fdm.h"

flux_face_FOU_vrans::flux_face_FOU_vrans(lexer *pp)
{
    p=pp;

}

flux_face_FOU_vrans::~flux_face_FOU_vrans()
{
}

void flux_face_FOU_vrans::u_flux(fdm* a,int ipol, field& uvel, double &uflux1, double &uflux2)
{
	if(ipol==1)
	{
        pip=1;
        
        if(p->flag1[Im1JK]>0)
        {
        if(0.5*(uvel(i,j,k)+uvel(i-1,j,k)) >= 0.0)
        uflux1 = uvel(i-1,j,k)*(1.0/a->porosity(i-1,j,k));
        
        if(0.5*(uvel(i,j,k)+uvel(i-1,j,k)) < 0.0)
        uflux1 = uvel(i,j,k)*(1.0/a->porosity(i-1,j,k));
        }
        
        if(p->flag1[Im1JK]<0)
        uflux1= 0.5*(uvel(i,j,k)+uvel(i-1,j,k))*(1.0/a->porosity(i-1,j,k));
        
        
        if(p->flag1[Ip1JK]>0)
        {
        if(0.5*(uvel(i,j,k)+uvel(i+1,j,k)) >= 0.0)
        uflux2 = uvel(i,j,k)*(1.0/a->porosity(i,j,k));
        
        if(0.5*(uvel(i,j,k)+uvel(i+1,j,k)) < 0.0)
        uflux2 = uvel(i+1,j,k)*(1.0/a->porosity(i,j,k));
        }
        
        if(p->flag1[Ip1JK]<0)
        uflux2= 0.5*(uvel(i,j,k)+uvel(i+1,j,k))*(1.0/a->porosity(i,j,k));

        pip=0;
	}

	if(ipol==2)
	{
	pip=1;
	uflux1= 0.5*(uvel(i-1,j,k)+uvel(i-1,j+1,k))*(1.0/(0.25*(a->porosity(i-1,j,k)+a->porosity(i-1,j+1,k)+a->porosity(i,j,k)+a->porosity(i,j+1,k))));
	uflux2= 0.5*(uvel(i,j,k)+uvel(i,j+1,k))*(1.0/(0.25*(a->porosity(i,j,k)+a->porosity(i,j+1,k)+a->porosity(i+1,j,k)+a->porosity(i+1,j+1,k))));
	pip=0;
	}

	if(ipol==3)
	{
	pip=1;
	uflux1= 0.5*(uvel(i-1,j,k)+uvel(i-1,j,k+1))*(1.0/(0.25*(a->porosity(i-1,j,k)+a->porosity(i-1,j,k+1)+a->porosity(i,j,k)+a->porosity(i,j,k+1))));
	uflux2= 0.5*(uvel(i,j,k)+uvel(i,j,k+1))*(1.0/(0.25*(a->porosity(i,j,k)+a->porosity(i,j,k+1)+a->porosity(i+1,j,k)+a->porosity(i+1,j,k+1))));
	pip=0;
	}

	if(ipol==4)
	{
	uflux1= uvel(i-1,j,k)*(1.0/(0.5*(a->porosity(i,j,k)+a->porosity(i-1,j,k))));
	uflux2= uvel(i,j,k)*(1.0/(0.5*(a->porosity(i+1,j,k)+a->porosity(i,j,k))));
	}

}

void flux_face_FOU_vrans::v_flux(fdm* a, int ipol, field& vvel, double &vflux1, double &vflux2)
{
	if(ipol==1)
	{
	pip=2;
	vflux1= 0.5*(vvel(i,j-1,k)+vvel(i+1,j-1,k))*(1.0/(0.25*(a->porosity(i,j-1,k)+a->porosity(i+1,j-1,k)+a->porosity(i,j,k)+a->porosity(i+1,j,k))));
	vflux2= 0.5*(vvel(i,j,k)+vvel(i+1,j,k))*(1.0/(0.25*(a->porosity(i,j,k)+a->porosity(i+1,j,k)+a->porosity(i,j+1,k)+a->porosity(i+1,j+1,k))));
	pip=0;
	}

	if(ipol==2)
	{
        pip=2;
        
        if(p->flag2[IJm1K]>0)
        {
        if(0.5*(vvel(i,j,k)+vvel(i,j-1,k)) >= 0.0)
        vflux1 = vvel(i,j-1,k)*(1.0/a->porosity(i,j-1,k));
        
        if(0.5*(vvel(i,j,k)+vvel(i,j-1,k)) < 0.0)
        vflux1 = vvel(i,j,k)*(1.0/a->porosity(i,j-1,k));
        }
        
        if(p->flag2[IJm1K]<0)
        vflux1= 0.5*(vvel(i,j,k)+vvel(i,j-1,k))*(1.0/a->porosity(i,j-1,k));
            
            
            
        if(p->flag2[IJp1K]>0)
        {
        if(0.5*(vvel(i,j,k)+vvel(i,j+1,k)) >= 0.0)
        vflux2 = vvel(i,j,k)*(1.0/a->porosity(i,j,k));
        
        if(0.5*(vvel(i,j,k)+vvel(i,j+1,k)) < 0.0)
        vflux2 = vvel(i,j+1,k)*(1.0/a->porosity(i,j,k));
        }
        
        if(p->flag2[IJp1K]<0)
        vflux2= 0.5*(vvel(i,j,k)+vvel(i,j+1,k))*(1.0/a->porosity(i,j,k));
        
        pip=0;
	}

	if(ipol==3)
	{
	pip=2;
	vflux1= 0.5*(vvel(i,j-1,k)+vvel(i,j-1,k+1))*(1.0/(0.25*(a->porosity(i,j-1,k)+a->porosity(i,j-1,k+1)+a->porosity(i,j,k)+a->porosity(i,j,k+1))));
	vflux2= 0.5*(vvel(i,j,k)+vvel(i,j,k+1))*(1.0/(0.25*(a->porosity(i,j,k)+a->porosity(i,j,k+1)+a->porosity(i,j+1,k)+a->porosity(i,j+1,k+1))));
	pip=0;
	}

	if(ipol==4)
	{
    pip=2;
	vflux1= vvel(i,j-1,k)*(1.0/(0.5*(a->porosity(i,j,k)+a->porosity(i,j-1,k))));
	vflux2= vvel(i,j,k)*(1.0/(0.5*(a->porosity(i,j+1,k)+a->porosity(i,j,k))));
    pip=0;
	}
}

void flux_face_FOU_vrans::w_flux(fdm* a, int ipol, field& wvel, double &wflux1, double &wflux2)
{

	if(ipol==1)
	{
	pip=3;
	wflux1= 0.5*(wvel(i,j,k-1)+wvel(i+1,j,k-1))*(1.0/(0.25*(a->porosity(i,j,k-1)+a->porosity(i+1,j,k-1)+a->porosity(i,j,k)+a->porosity(i+1,j,k))));
	wflux2= 0.5*(wvel(i,j,k)+wvel(i+1,j,k))*(1.0/(0.25*(a->porosity(i,j,k)+a->porosity(i+1,j,k)+a->porosity(i,j,k+1)+a->porosity(i+1,j,k+1))));
	pip=0;
	}

	if(ipol==2)
	{
	pip=3;
	wflux1= 0.5*(wvel(i,j,k-1)+wvel(i,j+1,k-1))*(1.0/(0.25*(a->porosity(i,j,k-1)+a->porosity(i,j+1,k-1)+a->porosity(i,j,k)+a->porosity(i,j+1,k))));
	wflux2= 0.5*(wvel(i,j,k)+wvel(i,j+1,k))*(1.0/(0.25*(a->porosity(i,j,k)+a->porosity(i,j+1,k)+a->porosity(i,j,k+1)+a->porosity(i,j+1,k+1))));
	pip=0;
	}


	if(ipol==3)
	{
        pip=3;
        
        if(p->flag3[IJKm1]>0)
        {
        if(0.5*(wvel(i,j,k)+wvel(i,j,k-1))>=0.0)
        wflux1 = wvel(i,j,k-1)*(1.0/a->porosity(i,j,k-1));
        
        if(0.5*(wvel(i,j,k)+wvel(i,j,k-1))<0.0)
        wflux1 = wvel(i,j,k)*(1.0/a->porosity(i,j,k-1));
        }
        
        if(p->flag3[IJKm1]<0)
        wflux1= 0.5*(wvel(i,j,k)+wvel(i,j,k-1))*(1.0/a->porosity(i,j,k-1));
        
        
        if(p->flag3[IJKp1]>0)
        {
        if(0.5*(wvel(i,j,k)+wvel(i,j,k+1))>=0.0)
        wflux2 = wvel(i,j,k)*(1.0/a->porosity(i,j,k));
        
        if(0.5*(wvel(i,j,k)+wvel(i,j,k))<0.0)
        wflux2 = wvel(i,j,k+1)*(1.0/a->porosity(i,j,k));
        }
        
        if(p->flag3[IJKp1]<0)
        wflux2= 0.5*(wvel(i,j,k)+wvel(i,j,k+1))*(1.0/a->porosity(i,j,k));

        pip=0;
	}

	if(ipol==4)
	{
    pip=3;
	wflux1= wvel(i,j,k-1)*(1.0/(0.5*(a->porosity(i,j,k)+a->porosity(i,j,k-1))));
	wflux2= wvel(i,j,k)*(1.0/(0.5*(a->porosity(i,j,k+1)+a->porosity(i,j,k))));
    pip=0;
	}
}
