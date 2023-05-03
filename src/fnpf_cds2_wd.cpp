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

#include"fnpf_cds2_wd.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"field.h"
#include"slice.h"
#include"vec.h"

fnpf_cds2_wd::fnpf_cds2_wd(lexer* p, fdm_fnpf *cc)
{
    c=cc;
}

fnpf_cds2_wd::~fnpf_cds2_wd()
{
}

double fnpf_cds2_wd::fx(lexer *p, field &f, double ivel1, double ivel2)
{
    if(p->wet[Im1J]>0 && p->wet[Ip1J]>0)
    return (f(i+1,j,k)-f(i-1,j,k))/(p->DXN[IP]+p->DXN[IM1]);
    
    else
    return 0.0;
}

double fnpf_cds2_wd::fy(lexer *p, field &f, double jvel1, double jvel2)
{
    if(p->wet[IJm1]>0 && p->wet[IJp1]>0)
    return (f(i,j+1,k)-f(i,j-1,k))/(p->DYN[JP]+p->DYN[JM1]);
    
    else
    return 0.0;
}

double fnpf_cds2_wd::fz(lexer *p, field &f, double kvel1, double kvel2)
{
    return (-1.5*f(i,j,k+1) + 2.0*f(i,j,k) - 0.5*f(i,j,k-1))/(-1.5*p->ZP[KP1] + 2.0*p->ZP[KP] - 0.5*p->ZP[KM1]);
}

double fnpf_cds2_wd::sx(lexer *p, slice &f, double ivel)
{
    if(p->wet[Im1J]>0 && p->wet[Ip1J]>0)
    return (f(i+1,j)-f(i-1,j))/(p->DXN[IP]+p->DXN[IM1]);
    
    else
    return 0.0;
}

double fnpf_cds2_wd::sy(lexer *p, slice &f, double jvel)
{
    if(p->wet[IJm1]>0 && p->wet[IJp1]>0)
    return (f(i,j+1)-f(i,j-1))/(p->DYN[JP]+p->DYN[JM1]);    
    
    else
    return 0.0;
}

double fnpf_cds2_wd::sz(lexer *p, double *f)
{
    return (-1.5*f[FIJK] + 2.0*f[FIJKm1] - 0.5*f[FIJKm2])/(-1.5*p->ZN[KP] + 2.0*p->ZN[KM1] - 0.5*p->ZN[KM2]);
}



