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

#include"flux_HJ_CDS4.h"
#include"lexer.h"
#include"fdm.h"

flux_HJ_CDS4::flux_HJ_CDS4(lexer *p)
{

}

flux_HJ_CDS4::~flux_HJ_CDS4()
{
}

void flux_HJ_CDS4::u_flux(fdm* a,int ipol, field& uvel, double &uflux1, double &uflux2)
{
	if(ipol==1)
	{
    pip=1;
	uflux1 = uvel(i,j,k);
    pip=0;
	}

	if(ipol==2)
	{
	pip=1;
	uflux1 = 0.25*(uvel(i,j,k) + uvel(i,j+1,k) + uvel(i-1,j,k) + uvel(i-1,j+1,k));
	pip=0;
	}

	if(ipol==3)
	{
	pip=1;
	uflux1 = 0.25*(uvel(i,j,k) + uvel(i,j,k+1) + uvel(i-1,j,k) + uvel(i-1,j,k+1));
	pip=0;
	}

	if(ipol==4)
	{
    pip=1;
	uflux1 = (1.0/16.0)*(-uvel(i-2,j,k)+9.0*uvel(i-1,j,k)+9.0*uvel(i,j,k)-uvel(i+1,j,k));
	pip=0;
	}
}

void flux_HJ_CDS4::v_flux(fdm* a,int ipol,field&vvel, double &vflux1, double &vflux2)
{
	if(ipol==1)
	{
	pip=2;
	vflux1 = 0.25*(vvel(i,j,k) + vvel(i+1,j,k) + vvel(i,j-1,k) + vvel(i+1,j-1,k));
	pip=0;
	}

	if(ipol==2)
	{
    pip=2;
	vflux1 = vvel(i,j,k);
    pip=0;
	}

	if(ipol==3)
	{
	pip=2;
	vflux1 = 0.25*(vvel(i,j,k) + vvel(i,j,k+1) + vvel(i,j-1,k) + vvel(i,j-1,k+1));
	pip=0;
	}

	if(ipol==4)
	{
    pip=2;
	vflux1 =  (1.0/16.0)*(-vvel(i,j-2,k)+9.0*vvel(i,j-1,k)+9.0*vvel(i,j,k)-vvel(i,j+1,k));
    pip=0;
	}
}

void flux_HJ_CDS4::w_flux(fdm* a,int ipol,field& wvel, double &wflux1, double &wflux2)
{

	if(ipol==1)
	{
	pip=3;
	wflux1 = 0.25*(wvel(i,j,k) + wvel(i+1,j,k) + wvel(i+1,j,k-1) + wvel(i,j,k-1));
	pip=0;
	}

	if(ipol==2)
	{
	pip=3;
	wflux1 = 0.25*(wvel(i,j,k) + wvel(i,j+1,k) + wvel(i,j+1,k-1) + wvel(i,j,k-1));
	pip=0;
	}

	if(ipol==3)
	{
    pip=3;
	wflux1 = wvel(i,j,k);
    pip=0;
	}

	if(ipol==4)
	{
    pip=3;
	wflux1 = (1.0/16.0)*(-wvel(i,j,k-2)+9.0*wvel(i,j,k-1)+9.0*wvel(i,j,k)-wvel(i,j,k+1));
    pip=0;
	}
}