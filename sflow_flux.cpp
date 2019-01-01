/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
/*
#include"sflow_flux.h"
#include"lexer.h"
#include"slice.h"

sflow_flux::sflow_flux(lexer *p)
{
}

sflow_flux::~sflow_flux()
{
}

void sflow_flux::iflux(int ipol, slice& uvel)
{
	if(ipol==1)
	{
    pip=1;
	ivel1= 0.5*(uvel(i,j)+uvel(i-1,j));
	ivel2= 0.5*(uvel(i,j)+uvel(i+1,j));
	pip=0;
	}

	if(ipol==2)
	{
	pip=1;
	ivel1= 0.5*(uvel(i-1,j)+uvel(i-1,j+1));
	ivel2= 0.5*(uvel(i,j)+uvel(i,j+1));
	pip=0;
	}

	if(ipol==4)
	{
	pip=1;
	ivel1= uvel(i-1,j);
	ivel2= uvel(i,j);
	pip=0;
	}
}

void sflow_flux::jflux(int ipol, slice& vvel)
{
	if(ipol==1)
	{
	pip=2;
	jvel1= 0.5*(vvel(i,j-1)+vvel(i+1,j-1));
	jvel2= 0.5*(vvel(i,j)+vvel(i+1,j));
	pip=0;
	}

	if(ipol==2)
	{
    pip=2;
	jvel1= 0.5*(vvel(i,j)+vvel(i,j-1));
	jvel2= 0.5*(vvel(i,j)+vvel(i,j+1));
	pip=0;
	}


	if(ipol==4)
	{
	pip=2;
	jvel1= vvel(i,j-1);
	jvel2= vvel(i,j);
	pip=0;
	}
}

void sflow_flux::ifluxC(fdm2D *b, int ipol, slice& uvel)
{
	if(ipol==1)
	{
    pip=1;
	ivel1= 0.25*(uvel(i,j) + uvel(i-1,j))*(b->hx(i,j) + b->hx(i-1,j));
	ivel2= 0.25*(uvel(i,j) + uvel(i+1,j))*(b->hx(i,j) + b->hx(i+1,j));
	pip=0;
	}

	if(ipol==2)
	{
	pip=1;
	ivel1= 0.25*(uvel(i-1,j) + uvel(i-1,j+1))*(b->hx(i-1,j) + b->hx(i-1,j+1));
	ivel2= 0.25*(uvel(i,j) + uvel(i,j+1))*(b->hx(i,j) + b->hx(i,j+1));
	pip=0;
	}

	if(ipol==4)
	{
	pip=1;
	ivel1= uvel(i-1,j)*b->hx(i-1,j);
	ivel2= uvel(i,j)*b->hx(i,j);
	pip=0;
	}
}

void sflow_flux::jfluxC(fdm2D *b, int ipol, slice& vvel)
{
	if(ipol==1)
	{
	pip=2;
	jvel1= 0.25*(vvel(i,j-1) + vvel(i+1,j-1))*(b->hy(i,j-1) + b->hy(i+1,j-1));
	jvel2= 0.25*(vvel(i,j) + vvel(i+1,j))*(b->hy(i,j) + b->hy(i+1,j));
	pip=0;
	}

	if(ipol==2)
	{
    pip=2;
	jvel1= 0.25*(vvel(i,j)+vvel(i,j-1))*(b->hy(i,j)+b->hy(i,j-1));
	jvel2= 0.25*(vvel(i,j)+vvel(i,j+1))*(b->hy(i,j)+b->hy(i,j+1));
	pip=0;
	}


	if(ipol==4)
	{
	pip=2;
	jvel1= vvel(i,j-1)*b->hy(i,j-1);
	jvel2= vvel(i,j)*b->hy(i,j);
	pip=0;
	}
}
 * 
 * 
 * 
 * 
 */ 

/*
void sflow_flux::ifluxC(fdm2D *b, int ipol, slice& uvel)
{
	if(ipol==1)
	{
    pip=1;
	ivel1= 0.5*(uvel(i,j)*b->hx(i,j) + uvel(i-1,j)*b->hx(i-1,j));
	ivel2= 0.5*(uvel(i,j)*b->hx(i,j) + uvel(i+1,j)*b->hx(i+1,j));
	pip=0;
	}

	if(ipol==2)
	{
	pip=1;
	ivel1= 0.5*(uvel(i-1,j)*b->hx(i-1,j) + uvel(i-1,j+1)*b->hx(i-1,j+1));
	ivel2= 0.5*(uvel(i,j)*b->hx(i,j) + uvel(i,j+1)*b->hx(i,j+1));
	pip=0;
	}

	if(ipol==4)
	{
	pip=1;
	ivel1= uvel(i-1,j);
	ivel2= uvel(i,j);
	pip=0;
	}
}

void sflow_flux::jfluxC(fdm2D *b, int ipol, slice& vvel)
{
	if(ipol==1)
	{
	pip=2;
	jvel1= 0.5*(vvel(i,j-1)*b->hy(i,j-1) + vvel(i+1,j-1)*b->hy(i+1,j-1));
	jvel2= 0.5*(vvel(i,j)*b->hy(i,j) + vvel(i+1,j)*b->hy(i+1,j));
	pip=0;
	}

	if(ipol==2)
	{
    pip=2;
	jvel1= 0.5*(vvel(i,j)*b->hy(i,j)+vvel(i,j-1)*b->hy(i,j-1));
	jvel2= 0.5*(vvel(i,j)*b->hy(i,j)+vvel(i,j+1)*b->hy(i,j+1));
	pip=0;
	}


	if(ipol==4)
	{
	pip=2;
	jvel1= vvel(i,j-1);
	jvel2= vvel(i,j);
	pip=0;
	}
}*/


/*
void sflow_flux::ifluxHJ(int ipol, slice &uvel)
{
	if(ipol==1)
	{
	iadvec= uvel(i,j);
	}

	if(ipol==2)
	{
	pip=1;
	iadvec=0.25*(uvel(i,j) + uvel(i,j+1) + uvel(i-1,j) + uvel(i-1,j+1));
	pip=0;
	}

	if(ipol==4)
	{
    pip=1;
	iadvec = 0.5*(uvel(i,j) + uvel(i-1,j));
	pip=0;
	}
}

void sflow_flux::jfluxHJ(int ipol, slice &vvel)
{
	if(ipol==1)
	{
	pip=2;
	jadvec=0.25*(vvel(i,j) + vvel(i+1,j) + vvel(i,j-1) + vvel(i+1,j-1));
	pip=0;
	}

	if(ipol==2)
	{
	jadvec= vvel(i,j);
	}

	if(ipol==4)
	{
    pip=2;
	jadvec = 0.5*(vvel(i,j) + vvel(i,j-1));
    pip=0;
	}
}*/
