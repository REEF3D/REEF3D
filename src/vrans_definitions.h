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

// NHFLOW VRANS conventions
// - U,V,W (and UH,VH,WH) are Darcy (filtered) velocities, pressure and eta are intrinsic.
// - Added mass follows van Gent (1995): (1 + c) du/dt = RHS,  c = gamma (1-n)/n^2  (gamma = B 260),
//   i.e. every momentum tendency, including the non-hydrostatic pressure correction, is scaled
//   by CPORNH = 1/(1+c).
//   NOTE: REEF3D::CFD uses CPOR = 1/(1+gamma) inside the porous zone (looping.h), a
//   porosity-independent added mass. The two modules therefore differ for the same B 260.
// - Continuity: n_s deta/dt + div(Q) = 0 with n_s the porosity of the surface cell.

#define PORVALNH  d->POR[IJK]

#define PORVALNH1  (0.5*(d->POR[IJK] + d->POR[Ip1JK]))
#define PORVALNH1m  (0.5*(d->POR[IJK] + d->POR[Im1JK]))
#define PORVALNH2  (0.5*(d->POR[IJK] + d->POR[IJp1K]))
#define PORVALNH2m  (0.5*(d->POR[IJK] + d->POR[IJm1K]))
#define PORVALNH3  (0.5*(d->POR[IJK] + d->POR[IJKp1]))

#define CPORNHval(n)  (1.0/(1.0 + (p->B260*(1.0-(n))/((n)*(n)))))

#define CPORNH  (1.0/(1.0 + (p->B260*(1.0-PORVALNH)/(PORVALNH*PORVALNH))))

#define CPORNH1m  (1.0/(1.0 + (p->B260*(1.0-PORVALNH1m)/(PORVALNH1m*PORVALNH1m))))
#define CPORNH1  (1.0/(1.0 + (p->B260*(1.0-PORVALNH1)/(PORVALNH1*PORVALNH1))))

#define CPORNH2m  (1.0/(1.0 + (p->B260*(1.0-PORVALNH2m)/(PORVALNH2m*PORVALNH2m))))
#define CPORNH2  (1.0/(1.0 + (p->B260*(1.0-PORVALNH2)/(PORVALNH2*PORVALNH2))))

/*
#define PORVALNH  1.0

#define PORVALNH1  1.0
#define PORVALNH1m  1.0
#define PORVALNH2  1.0
#define PORVALNH2m  1.0
#define PORVALNH3  1.0

#define CPORNHval(n)  1.0

#define CPORNH  1.0

#define CPORNH1m  1.0
#define CPORNH1  1.0

#define CPORNH2m  1.0
#define CPORNH2  1.0
*/
