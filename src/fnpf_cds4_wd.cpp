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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"fnpf_cds4_wd.h"
#include"lexer.h"
#include"field.h"
#include"slice.h"
#include"vec.h"

fnpf_cds4_wd::fnpf_cds4_wd(lexer* p)
{
}

fnpf_cds4_wd::~fnpf_cds4_wd()
{
}

double fnpf_cds4_wd::fx(lexer *p, field &f, double ivel1, double ivel2)
{
    return (-f(i+2,j,k) + 8.0*f(i+1,j,k) - 8.0*f(i-1,j,k) + f(i-2,j,k))
          /(-p->XP[IP2] + 8.0*p->XP[IP1] - 8.0*p->XP[IM1] + p->XP[IM2]);
}

double fnpf_cds4_wd::fy(lexer *p, field &f, double jvel1, double jvel2)
{
    return (-f(i,j+2,k) + 8.0*f(i,j+1,k) - 8.0*f(i,j-1,k) + f(i,j-2,k))
          /(-p->YP[JP2] + 8.0*p->YP[JP1] - 8.0*p->YP[JM1] + p->YP[JM2]);
}

double fnpf_cds4_wd::fz(lexer *p, field &f, double kvel1, double kvel2)
{
    return (-(25.0/12.0)*f(i,j,k+1) + 4.0*f(i,j,k) - 3.0*f(i,j,k-1) + (4.0/3.0)*f(i,j,k-2) - 0.25*f(i,j,k-3))
          /(-(25.0/12.0)*p->ZP[KP1] + 4.0*p->ZP[KP] - 3.0*p->ZP[KM1] + (4.0/3.0)*p->ZP[KM2] - 0.25*p->ZP[KM3]);
}

double fnpf_cds4_wd::sx(lexer *p, slice &f, double ivel)
{
    return (-f(i+2,j) + 8.0*f(i+1,j) - 8.0*f(i-1,j) + f(i-2,j))
          /(-p->XP[IP2] + 8.0*p->XP[IP1] - 8.0*p->XP[IM1] + p->XP[IM2]);
}

double fnpf_cds4_wd::sy(lexer *p, slice &f, double jvel)
{
    return (-f(i,j+2) + 8.0*f(i,j+1) - 8.0*f(i,j-1) + f(i,j-2))
          /(-p->YP[JP2] + 8.0*p->YP[JP1] - 8.0*p->YP[JM1] + p->YP[JM2]);   
}

double fnpf_cds4_wd::sz(lexer *p, double *f)
{
   // return (-(25.0/12.0)*f[FIJKp1] + 4.0*f[FIJK] - 3.0*f[FIJKm1] + (4.0/3.0)*f[FIJKm2] - 0.25*f[FIJKm3])
   //       /(-(25.0/12.0)*p->ZN[KP1] + 4.0*p->ZN[KP] - 3.0*p->ZN[KM1] + (4.0/3.0)*p->ZN[KM2] - 0.25*p->ZN[KM3]);
          
    return (-(25.0/12.0)*f[FIJK] + 4.0*f[FIJKm1] - 3.0*f[FIJKm2] + (4.0/3.0)*f[FIJKm3] - 0.25*f[FIJKm4])
          /(-(25.0/12.0)*p->ZN[KP] + 4.0*p->ZN[KM1] - 3.0*p->ZN[KM2] + (4.0/3.0)*p->ZN[KM3] - 0.25*p->ZN[KM4]);
          
        
}



