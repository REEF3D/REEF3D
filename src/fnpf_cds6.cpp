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

#include"fnpf_cds6.h"
#include"lexer.h"
#include"field.h"
#include"slice.h"
#include"vec.h"

fnpf_cds6::fnpf_cds6(lexer* p)
{
}

fnpf_cds6::~fnpf_cds6()
{
}

double fnpf_cds6::fx(lexer *p, field &f, double ivel1, double ivel2)
{
    return (f(i+3,j,k) - 9.0*f(i+2,j,k) + 45.0*f(i+1,j,k) - 45.0*f(i-1,j,k) + 9.0*f(i-2,j,k) - f(i-3,j,k))
          /(p->XP[IP3] - 9.0*p->XP[IP2] + 45.0*p->XP[IP1] - 45.0*p->XP[IM1] + 9.0*p->XP[IM2] - p->XP[IM3]);
}

double fnpf_cds6::fy(lexer *p, field &f, double jvel1, double jvel2)
{
    return (f(i,j+3,k) - 9.0*f(i,j+2,k) + 45.0*f(i,j+1,k) - 45.0*f(i,j-1,k) + 9.0*f(i,j-2,k) - f(i,j-3,k))
          /(p->YP[JP3] - 9.0*p->YP[JP3] + 45.0*p->YP[JP1] - 45.0*p->YP[JM1] + 9.0*p->YP[JM2] - p->YP[JM3]);
}

double fnpf_cds6::fz(lexer *p, field &f, double kvel1, double kvel2)
{
    return (-(49.0/20.0)*f(i,j,k+1) + 6.0*f(i,j,k) - 7.5*f(i,j,k-1) + (20.0/3.0)*f(i,j,k-2) - (15.0/4.0)*f(i,j,k-3) + (6.0/5.0)*f(i,j,k-4) - (1.0/6.0)*f(i,j,k-5))
          /(-(49.0/20.0)*p->ZP[KP1] + 6.0*p->ZP[KP] - 7.5*p->ZP[KM1] + (20.0/3.0)*p->ZP[KM2] - (15.0/4.0)*p->ZP[KM3] + (6.0/5.0)*p->ZP[KM4] - (1.0/6.0)*p->ZP[KM5]);
}

double fnpf_cds6::sx(lexer *p, slice &f, double ivel)
{
    return (f(i+3,j) - 9.0*f(i+2,j) + 45.0*f(i+1,j) - 45.0*f(i-1,j) + 9.0*f(i-2,j) - f(i-3,j))
          /(p->XP[IP3] - 9.0*p->XP[IP2] + 45.0*p->XP[IP1] - 45.0*p->XP[IM1] + 9.0*p->XP[IM2] - p->XP[IM3]);
}

double fnpf_cds6::sy(lexer *p, slice &f, double jvel)
{
    return (f(i,j+3) - 9.0*f(i,j+2) + 45.0*f(i,j+1) - 45.0*f(i,j-1) + 9.0*f(i,j-2) - f(i,j-3))
          /(p->YP[JP3] - 9.0*p->YP[JP3] + 45.0*p->YP[JP1] - 45.0*p->YP[JM1] + 9.0*p->YP[JM2] - p->YP[JM3]);  
}

double fnpf_cds6::sz(lexer *p, double *f)
{
    return (-(49.0/20.0)*f[FIJK] + 6.0*f[FIJKm1] - 7.5*f[FIJKm2] + (20.0/3.0)*f[FIJKm3] - 3.75*f[FIJKm4] + (6.0/5.0)*f[FIJKm5] - (1.0/6.0)*f[FIJKm6])
          /(-(49.0/20.0)*p->ZN[KP] + 6.0*p->ZN[KM1] - 7.5*p->ZN[KM2] + (20.0/3.0)*p->ZN[KM3] - 3.75*p->ZN[KM4] + (6.0/5.0)*p->ZN[KM5] - (1.0/6.0)*p->ZN[KM6]);
     /*     
    return (-(25.0/12.0)*f[FIJK] + 4.0*f[FIJKm1] - 3.0*f[FIJKm2] + (4.0/3.0)*f[FIJKm3] - 0.25*f[FIJKm4])
          /(-(25.0/12.0)*p->ZN[KP] + 4.0*p->ZN[KM1] - 3.0*p->ZN[KM2] + (4.0/3.0)*p->ZN[KM3] - 0.25*p->ZN[KM4]);*/
}



