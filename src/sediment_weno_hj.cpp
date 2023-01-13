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

#include"sediment_weno_hj.h"
#include"lexer.h"
#include"vec.h"
#include"fnpf_discrete_weights.h"

sediment_weno_hj::sediment_weno_hj(lexer* p) :  ddweno_f_nug(p)
{
    p->Darray(ckz,p->knoz+1+4*marge,5);
    
    fnpf_discrete_weights dw(p);

    dw.ck_weights(p, ckz, p->ZN, p->knoz+1, 1, 4, 6);
}

sediment_weno_hj::~sediment_weno_hj()
{
}
double sediment_weno_hj::sx(lexer *p, slice &f, double ivel1, double ivel2)
{
    grad=0.0;
    
    ivel = 0.5*(ivel1+ivel2);
    
    if(ivel>0.0)
    grad=dswenox(f,1.0)*ivel;
    
    if(ivel<0.0)
    grad=dswenox(f,-1.0)*ivel;
    
    return grad;
}

double sediment_weno_hj::sy(lexer *p, slice &f, double jvel1, double jvel2)
{
    grad=0.0;
    
    jvel = 0.5*(jvel1+jvel2);
    
    if(jvel>0.0)
    grad=dswenoy(f,1.0)*jvel;
    
    if(jvel<0.0)
    grad=dswenoy(f,-1.0)*jvel;
    
    return grad;   
}


