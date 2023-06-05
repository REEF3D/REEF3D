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

#include"nhflow_flux_face_cds2.h"
#include"lexer.h"
#include"fdm.h"

nhflow_flux_face_cds2::nhflow_flux_face_cds2(lexer *pp)
{
p=pp;
}

nhflow_flux_face_cds2::~nhflow_flux_face_cds2()
{
}

void nhflow_flux_face_cds2::u_flux(fdm_nhf *d,int ipol, double *U, double &uflux1, double &uflux2)
{
    if(ipol==1)
	{
	uflux1= 0.5*(U[IJK]+U[Im1JK]);
	uflux2= 0.5*(U[IJK]+U[Ip1JK]);
	}

	if(ipol==2)
	{
	uflux1= 0.5*(U[Im1JK]+U[Im1Jp1K]);
	uflux2= 0.5*(U[IJK]+U[IJp1K]);
	}

	if(ipol==3)
	{
	uflux1= 0.5*(U[Im1JK]+U[Im1JKp1]);
	uflux2= 0.5*(U[IJK]+U[IJKp1]);
	}

	if(ipol==4)
	{
	uflux1= U[Im1JK];
	uflux2= U[IJK];
	}
}

void nhflow_flux_face_cds2::v_flux(fdm_nhf *d, int ipol, double *V, double &vflux1, double &vflux2)
{
    if(ipol==1)
	{
	vflux1= 0.5*(V[IJm1K]+V[Ip1Jm1K]);
	vflux2= 0.5*(V[IJK]+V[Ip1JK]);
	}

	if(ipol==2)
	{
	vflux1= 0.5*(V[IJK]+V[IJm1K]);
	vflux2= 0.5*(V[IJK]+V[IJp1K]);
	}

	if(ipol==3)
	{
	vflux1= 0.5*(V[IJp1K]+V[IJm1Kp1]);
	vflux2= 0.5*(V[IJK]+V[IJKp1]);
	}

	if(ipol==4)
	{
	vflux1= V[IJm1K];
	vflux2= V[IJK];
	}
}

void nhflow_flux_face_cds2::w_flux(fdm_nhf *d, int ipol, double *W, double &wflux1, double &wflux2)
{
    if(ipol==1)
	{
	pip=3;
	wflux1= 0.5*(W[FIJKm1]+W[FIp1JKm1]);
	wflux2= 0.5*(W[FIJK]+W[FIp1JK]);
	pip=0;
	}

	if(ipol==2)
	{
	pip=3;
	wflux1= 0.5*(W[FIJKm1]+W[FIJp1Km1]);
	wflux2= 0.5*(W[FIJK]+W[FIJp1K]);
	pip=0;
	}

	if(ipol==3)
	{
    pip=3;
	wflux1= 0.5*(W[FIJK]+W[FIJKm1]);
	wflux2= 0.5*(W[FIJK]+W[FIJKp1]);
	pip=0;
	}

	if(ipol==4)
	{
    pip=3;
	wflux1= W[FIJKm1];
	wflux2= W[FIJK];
    pip=0;
	}
    
    if(p->A501==2)
    {
    wflux1= W[FIJK];
	wflux2= W[FIJKp1];
    
    }
}

void nhflow_flux_face_cds2::omega_flux(lexer *p, fdm_nhf *d, int ipol, double *U, double *V, double *W, double &wflux1, double &wflux2)
{
    wflux1 =  p->sigt[FIJK]
            + 0.5*(U[IJK]+U[IJKm1])
            + 0.5*(V[IJK]+V[IJKm1])
            + 0.5*(W[FIJK]+W[IJKm1])*p->sigz[IJ];
    
	wflux2 =  p->sigt[FIJKp1]
            + 0.5*(U[IJK]+U[IJKp1])
            + 0.5*(V[IJK]+V[IJKp1])
            + 0.5*(W[FIJK]+W[IJKp1])*p->sigz[IJ];
}

