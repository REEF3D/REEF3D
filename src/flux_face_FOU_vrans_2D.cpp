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

#include"flux_face_FOU_vrans_2D.h"
#include"lexer.h"
#include"fdm.h"

flux_face_FOU_vrans_2D::flux_face_FOU_vrans_2D(lexer *pp)
{
    p=pp;

}

flux_face_FOU_vrans_2D::~flux_face_FOU_vrans_2D()
{
}

void flux_face_FOU_vrans_2D::u_flux(fdm* a,int ipol, field& uvel, double &uflux1, double &uflux2)
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
	uflux1= uvel(i-1,j,k)*(1.0/(0.5*(a->porosity(i-1,j,k)+a->porosity(i,j,k))));
	uflux2= uvel(i,j,k)*(1.0/(0.5*(a->porosity(i,j,k)+a->porosity(i+1,j,k))));
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

void flux_face_FOU_vrans_2D::v_flux(fdm* a, int ipol, field& vvel, double &vflux1, double &vflux2)
{
	vflux1= 0.0;
	vflux2= 0.0;
}

void flux_face_FOU_vrans_2D::w_flux(fdm* a, int ipol, field& wvel, double &wflux1, double &wflux2)
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
	wflux1= wvel(i,j,k-1)*(1.0/(0.5*(a->porosity(i,j,k-1)+a->porosity(i,j,k))));
	wflux2= wvel(i,j,k)*(1.0/(0.5*(a->porosity(i,j,k)+a->porosity(i,j,k+1))));
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
