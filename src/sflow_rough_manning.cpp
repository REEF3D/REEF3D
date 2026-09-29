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

#include"sflow_rough_manning.h"
#include"lexer.h"
#include"fdm2D.h"

#define WLIJ (WL(i,j)>p->A244?WL(i,j):1.0e20)

sflow_rough_manning::sflow_rough_manning(lexer* p) 
{
}

sflow_rough_manning::~sflow_rough_manning()
{
}

// Manning friction, tau_b/rho = g n^2 |u| u / h^(1/3), n = ks^(1/6)/20

void sflow_rough_manning::u_source(lexer *p, fdm2D *b, slice &U, slice &V, slice &WL)
{
    SLICELOOP4
    WETDRY
    {
    manning = pow(b->ks(i,j),1.0/6.0)/20.0;
    cf = pow(manning,2.0)*fabs(p->W22)/pow(WLIJ,1.0/3.0);
    
    b->F(i,j) -= cf*U(i,j)*sqrt(U(i,j)*U(i,j) + V(i,j)*V(i,j));
    }
}

void sflow_rough_manning::v_source(lexer *p, fdm2D *b, slice &U, slice &V, slice &WL)
{
    SLICELOOP4
    WETDRY
    {
    manning = pow(b->ks(i,j),1.0/6.0)/20.0;
    cf = pow(manning,2.0)*fabs(p->W22)/pow(WLIJ,1.0/3.0);
    
    b->G(i,j) -= cf*V(i,j)*sqrt(U(i,j)*U(i,j) + V(i,j)*V(i,j));
    }
}
