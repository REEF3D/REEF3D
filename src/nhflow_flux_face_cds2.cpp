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
	uflux1= 0.5*(U[IJK]+U[Im1JK]);
	uflux2= 0.5*(U[IJK]+U[Ip1JK]);
}

void nhflow_flux_face_cds2::v_flux(fdm_nhf *d, int ipol, double *V, double &vflux1, double &vflux2)
{
	vflux1= 0.5*(V[IJK]+V[IJm1K]);
	vflux2= 0.5*(V[IJK]+V[IJp1K]);
}

void nhflow_flux_face_cds2::w_flux(fdm_nhf *d, int ipol, double *W, double &wflux1, double &wflux2)
{
	wflux1= 0.5*(W[IJK]+W[IJKm1]);
	wflux2= 0.5*(W[IJK]+W[IJKp1]);
}

void nhflow_flux_face_cds2::omega_flux(lexer *p, fdm_nhf *d, int ipol, double *U, double *V, double *W, double &wflux1, double &wflux2)
{
    /*{
    a->omega[IJK] =  p->sigt[FIJKp1]
                    +  0.25*(U[Im1JK] + U[Im1JKp1] + U[IJK] + U[IJKp1])*p->sigx[FIJKp1]
                    +  0.25*(V[IJm1K] + V[IJm1Kp1] + V[IJK] + V[IJKp1])*p->sigy[FIJKp1]
                    +  W[IJK]*p->sigz[IJ];
    }*/
    
    wflux1 =  p->sigt[FIJK]
            + 0.5*(U[IJK]+U[IJKm1])
            + 0.5*(V[IJK]+V[IJKm1])
            + 0.5*(W[IJK]+W[IJKm1])*p->sigz[IJ];
    
	wflux2 =  p->sigt[FIJKp1]
            + 0.5*(U[IJK]+U[IJKp1])
            + 0.5*(V[IJK]+V[IJKp1])
            + 0.5*(W[IJK]+W[IJKp1])*p->sigz[IJ];
}

