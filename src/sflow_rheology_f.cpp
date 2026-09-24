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

#include"sflow_rheology_f.h"
#include"lexer.h"
#include"fdm2D.h"

sflow_rheology_f::sflow_rheology_f(lexer* p) 
{
    tanphi=0.0;
    if(p->W101>0)
    tanphi=tan(p->W102_phi*(PI/180.0));
}

sflow_rheology_f::~sflow_rheology_f()
{
}

void sflow_rheology_f::u_source(lexer *p, fdm2D *b, slice &U, slice &V, slice &WL)
{
    SLICELOOP4
    WETDRY
    {   
    u_abs = sqrt(U(i,j)*U(i,j) + V(i,j)*V(i,j));
    
    press = b->press(i,j) + p->W1*fabs(p->W22)*WL(i,j);
    
    tau_zx = bingham(p,b,U(i,j),u_abs,press,HPIJ);
    
    b->F(i,j) -= tau_zx/p->W1;
    }
}

void sflow_rheology_f::v_source(lexer *p, fdm2D *b, slice &U, slice &V, slice &WL)
{
    SLICELOOP4
    WETDRY
    {    
    u_abs = sqrt(U(i,j)*U(i,j) + V(i,j)*V(i,j));
    
    press = b->press(i,j) + p->W1*fabs(p->W22)*WL(i,j);
    
    tau_zy = bingham(p,b,V(i,j),u_abs,press,HPIJ);
    
    b->G(i,j) -= tau_zy/p->W1;
    }   
}
