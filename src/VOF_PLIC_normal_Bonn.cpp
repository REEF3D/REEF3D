/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
Author: Fabian Knoblauch
--------------------------------------------------------------------*/

#include"VOF_PLIC.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"solver.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"fluid_update_vof.h"
#include"heat.h"
#include"hires.h"
#include"weno_hj.h"

void VOF_PLIC::simpleNormal_Bonn
(
    fdm* a,
    lexer* p
)
{
    double nx_simp,ny_simp,nz_simp,nsum;
    nx_simp=(phistep(i+1,j,k)-phistep(i-1,j,k))/(p->DXP[IM1]+p->DXP[IP]);
    ny_simp=(phistep(i,j+1,k)-phistep(i,j-1,k))/(p->DYP[IM1]+p->DYP[IP]);
    nz_simp=(phistep(i,j,k+1)-phistep(i,j,k-1))/(p->DZP[IM1]+p->DZP[IP]);
    nsum=sqrt(nx_simp*nx_simp+ny_simp*ny_simp+nz_simp*nz_simp);
    nx(i,j,k)=nx_simp/nsum;
    ny(i,j,k)=ny_simp/nsum;
    nz(i,j,k)=nz_simp/nsum;
}