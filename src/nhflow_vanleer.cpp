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

#include"nhflow_vanleer.h"
#include"lexer.h"
#include"fdm.h"

nhflow_vanleer::nhflow_vanleer (lexer *pp)
{
    p=pp;
}

nhflow_vanleer::~nhflow_vanleer()
{
}

double nhflow_vanleer::iphi(double *F,int n1, int n2, int q1, int q2)
{
    denom=(F[(i+q1-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin] - F[(i+q2-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]);
    r=(F[(i+n1-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin] - F[(i+n2-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin])/(denom+1.0e-20);

    phi = (r+fabs(r))/(1.0+fabs(r));
    
    return phi;
}

double nhflow_vanleer::jphi(double *F,int n1, int n2, int q1, int q2)
{
    denom=(F[(i-p->imin)*p->jmax*p->kmax + (j+q1-p->jmin)*p->kmax + k-p->kmin] - F[(i-p->imin)*p->jmax*p->kmax + (j+q2-p->jmin)*p->kmax + k-p->kmin]);
    r=(F[(i-p->imin)*p->jmax*p->kmax + (j+n1-p->jmin)*p->kmax + k-p->kmin] - F[(i-p->imin)*p->jmax*p->kmax + (j+n2-p->jmin)*p->kmax + k-p->kmin])/(denom+1.0e-20);

    phi = (r+fabs(r))/(1.0+fabs(r));
	
    return phi;
}

double nhflow_vanleer::kphi(double *F,int n1, int n2, int q1, int q2)
{
    denom=(F[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k+q1-p->kmin] - F[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k+q2-p->kmin]);
    r=(F[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k+n1-p->kmin] - F[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k+n2-p->kmin])/(denom+1.0e-20);

    phi = (r+fabs(r))/(1.0+fabs(r));
	
    return phi;
}
