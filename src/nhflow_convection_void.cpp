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

#include"nhflow_convection_void.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"nhflow_flux_face_cds2.h"

nhflow_convection_void::nhflow_convection_void(lexer* p)
{
}

nhflow_convection_void::~nhflow_convection_void()
{
}

void nhflow_convection_void::precalc(lexer* p, fdm_nhf* d, double *F, int ipol, double *U, double *V, double *W, slice &eta)
{
}

void nhflow_convection_void::start(lexer* p, fdm_nhf* d, double *B, int ipol, double *U, double *V, double *W, slice &eta)
{
    
}