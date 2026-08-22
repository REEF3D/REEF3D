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

#include"nhflow_scalar_advec_CDS2.h"
#include"lexer.h"

nhflow_scalar_advec_CDS2::nhflow_scalar_advec_CDS2(lexer *pp)
{
p=pp;
}

nhflow_scalar_advec_CDS2::~nhflow_scalar_advec_CDS2()
{
}

void nhflow_scalar_advec_CDS2::uadvec(int ipol, double *U, double &uflux1, double &uflux2)
{
	uflux1 = 0.5*(U[IJK]+U[Im1JK]);
	uflux2 = 0.5*(U[IJK]+U[Ip1JK]);
}

void nhflow_scalar_advec_CDS2::vadvec(int ipol, double *V, double &vflux1, double &vflux2)
{
	vflux1 = 0.5*(V[IJK]+V[IJm1K]);
	vflux2 = 0.5*(V[IJK]+V[IJp1K]);
}

void nhflow_scalar_advec_CDS2::wadvec(int ipol, double *W, double &wflux1, double &wflux2)
{
	wflux1 = 0.5*(W[IJK]+W[IJKm1]);
	wflux2 = 0.5*(W[IJK]+W[IJKp1]);
}


