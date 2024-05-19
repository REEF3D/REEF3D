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

void sflow_rheology_f::u_source(lexer *p, fdm2D *b, slice &u, slice &v)
{
    SLICELOOP1
    {   
    u_abs = sqrt(pow(u(i,j),2.0) + pow(0.25*(v(i,j-1) + v(i,j) + v(i+1,j-1) + v(i+1,j)),2.0));
    press = 0.5*(b->press(i,j) + b->press(i+1,j)) + p->W1*fabs(p->W22)*0.5*(b->hp(i,j) + b->hp(i+1,j));
    
    tau_zx = bingham(p,b,u(i,j),u_abs,press,HXIJ);
    
    
    
    b->F(i,j) -= tau_zx/(HXIJ*p->W1);
    }
}

void sflow_rheology_f::v_source(lexer *p, fdm2D *b, slice &u, slice &v)
{
    SLICELOOP2
    {    
    u_abs = sqrt(pow(0.25*(u(i-1,j) + u(i,j) + u(i-1,j+1) + u(i,j+1)),2.0) + pow(v(i,j),2.0));
    press = 0.5*(b->press(i,j) + b->press(i,j+1)) + p->W1*fabs(p->W22)*0.5*(b->hp(i,j) + b->hp(i,j+1));
    
    tau_zy = bingham(p,b,v(i,j),u_abs,press,HYIJ);
    
    b->G(i,j) -= tau_zy/(HYIJ*p->W1);
    }   
}

