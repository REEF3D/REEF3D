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

#include"sediment_cds.h"
#include"lexer.h"
#include"slice.h"
#include"vec.h"
#include"fnpf_discrete_weights.h"

sediment_cds::sediment_cds(lexer* p) 
{

}

sediment_cds::~sediment_cds()
{
}

double sediment_cds::sx(lexer *p, slice &f, double ivel1, double ivel2)
{   
    grad = (f(i+1,j)*ivel2-f(i-1,j)*ivel1)/(p->DXP[IP]+p->DXP[IM1]);
        
    return grad;
}

double sediment_cds::sy(lexer *p, slice &f, double jvel1, double jvel2)
{
    grad = (f(i,j+1)*jvel2-f(i,j-1)*jvel1)/(p->DYP[JP]+p->DYP[JM1]);
			  
    return grad;  
}

