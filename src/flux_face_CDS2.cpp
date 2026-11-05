/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"flux_face_CDS2.h"
#include"lexer.h"
#include"fdm.h"

flux_face_CDS2::flux_face_CDS2(lexer *pp)
{
p=pp;
}

flux_face_CDS2::~flux_face_CDS2()
{
}

void flux_face_CDS2::u_flux(fdm* a,int ipol, field& uvel, double &uflux1, double &uflux2)
{
	if(ipol==1)
	{
	uflux1= 0.5*(uvel(i,j,k)+uvel(i-1,j,k));
	uflux2= 0.5*(uvel(i,j,k)+uvel(i+1,j,k));
	}

	if(ipol==2)
	{
	uflux1= 0.5*(uvel(i-1,j,k)+uvel(i-1,j+1,k));
	uflux2= 0.5*(uvel(i,j,k)+uvel(i,j+1,k));
	}

	if(ipol==3)
	{
	uflux1= 0.5*(uvel(i-1,j,k)+uvel(i-1,j,k+1));
	uflux2= 0.5*(uvel(i,j,k)+uvel(i,j,k+1));
	}

	if(ipol==4)
	{
	uflux1= uvel(i-1,j,k);
	uflux2= uvel(i,j,k);
	}
}

void flux_face_CDS2::v_flux(fdm* a, int ipol, field& vvel, double &vflux1, double &vflux2)
{
	if(ipol==1)
	{
	vflux1= 0.5*(vvel(i,j-1,k)+vvel(i+1,j-1,k));
	vflux2= 0.5*(vvel(i,j,k)+vvel(i+1,j,k));
	}

	if(ipol==2)
	{
	vflux1= 0.5*(vvel(i,j,k)+vvel(i,j-1,k));
	vflux2= 0.5*(vvel(i,j,k)+vvel(i,j+1,k));
	}

	if(ipol==3)
	{
	vflux1= 0.5*(vvel(i,j-1,k)+vvel(i,j-1,k+1));
	vflux2= 0.5*(vvel(i,j,k)+vvel(i,j,k+1));
	}

	if(ipol==4)
	{
	vflux1= vvel(i,j-1,k);
	vflux2= vvel(i,j,k);
	}
}

void flux_face_CDS2::w_flux(fdm* a, int ipol, field& wvel, double &wflux1, double &wflux2)
{

	if(ipol==1)
	{
	wflux1= 0.5*(wvel(i,j,k-1)+wvel(i+1,j,k-1));
	wflux2= 0.5*(wvel(i,j,k)+wvel(i+1,j,k));
	}

	if(ipol==2)
	{
	wflux1= 0.5*(wvel(i,j,k-1)+wvel(i,j+1,k-1));
	wflux2= 0.5*(wvel(i,j,k)+wvel(i,j+1,k));
	}

	if(ipol==3)
	{
	wflux1= 0.5*(wvel(i,j,k)+wvel(i,j,k-1));
	wflux2= 0.5*(wvel(i,j,k)+wvel(i,j,k+1));
	}

	if(ipol==4)
	{
	wflux1= wvel(i,j,k-1);
	wflux2= wvel(i,j,k);
	}
    
    // no vertical flux next to a solid or the topography (static geometry).  Not next to a
    // floating body (flagsf4 < 0 in fb < 0 too): the x and y fluxes are not cut there, and a body
    // moving in z (a propeller turning about x) then builds up a spurious jet in w at its surface
    if(p->flagsf4[IJKm1]<0 && static_geometry(a,i,j,k-1))
    wflux1=0.0;
    
    if(p->flagsf4[IJKp1]<0 && static_geometry(a,i,j,k+1))
    wflux2=0.0;
}

bool flux_face_CDS2::static_geometry(fdm *a, int ii, int jj, int kk)
{
    if(p->X10!=1)
    return true;
    
    return (p->solidread==1 && a->solid(ii,jj,kk)<0.0) || (p->toporead==1 && a->topo(ii,jj,kk)<0.0);
}


