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

#include"flux_face_CDS2_vrans_2D.h"
#include"lexer.h"
#include"fdm.h"

flux_face_CDS2_vrans_2D::flux_face_CDS2_vrans_2D(lexer *p)
{

}

flux_face_CDS2_vrans_2D::~flux_face_CDS2_vrans_2D()
{
}

void flux_face_CDS2_vrans_2D::u_flux(fdm* a,int ipol, field& uvel, double &uflux1, double &uflux2)
{
	if(ipol==1)
	{
	uflux1= 0.5*(uvel(i,j,k)+uvel(i-1,j,k))*(1.0/a->porosity(i-1,j,k));
	uflux2= 0.5*(uvel(i,j,k)+uvel(i+1,j,k))*(1.0/a->porosity(i,j,k));
	}

	if(ipol==2)
	{
	uflux1= uvel(i-1,j,k)*(1.0/(0.5*(a->porosity(i-1,j,k)+a->porosity(i,j,k))));
	uflux2= uvel(i,j,k)*(1.0/(0.5*(a->porosity(i,j,k)+a->porosity(i+1,j,k))));
	}

	if(ipol==3)
	{
	uflux1= 0.5*(uvel(i-1,j,k)+uvel(i-1,j,k+1))*(1.0/(0.25*(a->porosity(i-1,j,k)+a->porosity(i-1,j,k+1)+a->porosity(i,j,k)+a->porosity(i,j,k+1))));
	uflux2= 0.5*(uvel(i,j,k)+uvel(i,j,k+1))*(1.0/(0.25*(a->porosity(i,j,k)+a->porosity(i,j,k+1)+a->porosity(i+1,j,k)+a->porosity(i+1,j,k+1))));
	}

	if(ipol==4)
	{
	uflux1= uvel(i-1,j,k)*(1.0/(0.5*(a->porosity(i,j,k)+a->porosity(i-1,j,k))));
	uflux2= uvel(i,j,k)*(1.0/(0.5*(a->porosity(i+1,j,k)+a->porosity(i,j,k))));
	}

}

void flux_face_CDS2_vrans_2D::v_flux(fdm* a, int ipol, field& vvel, double &vflux1, double &vflux2)
{
	vflux1= 0.0;
	vflux2= 0.0;
}

void flux_face_CDS2_vrans_2D::w_flux(fdm* a, int ipol, field& wvel, double &wflux1, double &wflux2)
{

	if(ipol==1)
	{
	wflux1= 0.5*(wvel(i,j,k-1)+wvel(i+1,j,k-1))*(1.0/(0.25*(a->porosity(i,j,k-1)+a->porosity(i+1,j,k-1)+a->porosity(i,j,k)+a->porosity(i+1,j,k))));
	wflux2= 0.5*(wvel(i,j,k)+wvel(i+1,j,k))*(1.0/(0.25*(a->porosity(i,j,k)+a->porosity(i+1,j,k)+a->porosity(i,j,k+1)+a->porosity(i+1,j,k+1))));
	}

	if(ipol==2)
	{
	wflux1= wvel(i,j,k-1)*(1.0/(0.5*(a->porosity(i,j,k-1)+a->porosity(i,j,k))));
	wflux2= wvel(i,j,k)*(1.0/(0.5*(a->porosity(i,j,k)+a->porosity(i,j,k+1))));
	}

	if(ipol==3)
	{
	wflux1= 0.5*(wvel(i,j,k)+wvel(i,j,k-1))*(1.0/a->porosity(i,j,k-1));
	wflux2= 0.5*(wvel(i,j,k)+wvel(i,j,k+1))*(1.0/a->porosity(i,j,k));
	}

	if(ipol==4)
	{
	wflux1= wvel(i,j,k-1)*(1.0/(0.5*(a->porosity(i,j,k)+a->porosity(i,j,k-1))));
	wflux2= wvel(i,j,k)*(1.0/(0.5*(a->porosity(i,j,k+1)+a->porosity(i,j,k))));
	}
}
