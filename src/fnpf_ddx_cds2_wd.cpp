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

#include"fnpf_ddx_cds2_wd.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"slice.h"

fnpf_ddx_cds2_wd::fnpf_ddx_cds2_wd(lexer* p, fdm_fnpf *cc)
{
    c=cc;
}

fnpf_ddx_cds2_wd::~fnpf_ddx_cds2_wd()
{
}

double fnpf_ddx_cds2_wd::sxx(lexer *p, slice &f)
{
    if(p->wet[Im1J]>0 && p->wet[Ip1J]>0)
    return ((f(i+1,j)-f(i,j))/p->DXP[IP] - (f(i,j)-f(i-1,j))/p->DXP[IM1])/p->DXN[IP];
    
    else
    return 0.0;
}

double fnpf_ddx_cds2_wd::syy(lexer *p, slice &f)
{
    if(p->wet[IJm1]>0 && p->wet[IJp1]>0)
    return ((f(i,j+1)-f(i,j))/p->DYP[JP] - (f(i,j)-f(i,j-1))/p->DYP[JM1])/p->DYN[JP];   

    else
    return 0.0; 
}



