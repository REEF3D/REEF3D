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

#include"nhflow_hxy_flux.h"
#include"lexer.h"
#include"slice.h"

nhflow_hxy_flux::nhflow_hxy_flux(lexer *p)
{
}

nhflow_hxy_flux::~nhflow_hxy_flux()
{
}

void nhflow_hxy_flux::u_flux(int ipol, slice& uvel, double &uflux1, double &uflux2)
{
	if(ipol==1)
	{
    pip=1;
	uflux1 = 0.5*(uvel(i,j)+uvel(i-1,j));
	uflux2 = 0.5*(uvel(i,j)+uvel(i+1,j));
	pip=0;
	}

	if(ipol==2)
	{
	pip=1;
	uflux1 = 0.5*(uvel(i-1,j)+uvel(i-1,j+1));
	uflux2 = 0.5*(uvel(i,j)+uvel(i,j+1));
	pip=0;
	}

	if(ipol==4)
	{
	pip=1;
	uflux1 = uvel(i-1,j);
	uflux2 = uvel(i,j);
	pip=0;
	}
}

void nhflow_hxy_flux::v_flux(int ipol, slice& vvel, double &vflux1, double &vflux2)
{
	if(ipol==1)
	{
	pip=2;
	vflux1 = 0.5*(vvel(i,j-1)+vvel(i+1,j-1));
	vflux2 = 0.5*(vvel(i,j)+vvel(i+1,j));
	pip=0;
	}

	if(ipol==2)
	{
    pip=2;
	vflux1 = 0.5*(vvel(i,j)+vvel(i,j-1));
	vflux2 = 0.5*(vvel(i,j)+vvel(i,j+1));
	pip=0;
	}


	if(ipol==4)
	{
	pip=2;
	vflux1 = vvel(i,j-1);
	vflux2 = vvel(i,j);
	pip=0;
	}
}
