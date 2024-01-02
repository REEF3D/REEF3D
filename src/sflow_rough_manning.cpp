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

#include"sflow_rough_manning.h"
#include"lexer.h"
#include"fdm2D.h"

#define HXIJ (fabs(b->hx(i,j))>1.0e-20?b->hx(i,j):1.0e20)
#define HYIJ (fabs(b->hy(i,j))>1.0e-20?b->hy(i,j):1.0e20)

sflow_rough_manning::sflow_rough_manning(lexer* p) 
{
}

sflow_rough_manning::~sflow_rough_manning()
{
}

void sflow_rough_manning::u_source(lexer *p, fdm2D *b, slice &u)
{
    SLICELOOP1
    {
    manning = pow(0.5*(b->ks(i,j)+b->ks(i+1,j)),1.0/6.0)/20.0;
    
    cf = pow(manning,2.0)*9.81/pow(HXIJ,1.0/3.0);
    
    b->F(i,j) -= cf*u(i,j)*fabs(u(i,j))*(1.0/HXIJ);
    }
}

void sflow_rough_manning::v_source(lexer *p, fdm2D *b, slice &v)
{
    SLICELOOP2
    {
    manning = pow(0.5*(b->ks(i,j)+b->ks(i,j+1)),1.0/6.0)/20.0;
    
    cf = pow(manning,2.0)*9.81/pow(HYIJ,1.0/3.0);
    
    b->G(i,j) -= cf*v(i,j)*fabs(v(i,j))*(1.0/HYIJ);
    }   
}

