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
    if(p->S31==1)
    grad = (0.5*(f(i,j)+f(i+1,j))*ivel2 - 0.5*(f(i-1,j)+f(i,j))*ivel1)/p->DXN[IP];
    
    // S 31 2/3: face form (conservative on non-uniform grids, equals the central difference on a uniform grid),
    // closed faces (zero face direction) carry no flux
    if(p->S31>=2)
    grad = (0.5*(f(i,j)+f(i+1,j))*(ivel2!=0.0?1.0:0.0) - 0.5*(f(i-1,j)+f(i,j))*(ivel1!=0.0?1.0:0.0))/p->DXN[IP];
        
    return grad;
}

double sediment_cds::sy(lexer *p, slice &f, double jvel1, double jvel2)
{
    if(p->S31==1)
    grad = (0.5*(f(i,j)+f(i,j+1))*jvel2 - 0.5*(f(i,j-1)+f(i,j))*jvel1)/p->DYN[JP];
    
    if(p->S31>=2)
    grad = (0.5*(f(i,j)+f(i,j+1))*(jvel2!=0.0?1.0:0.0) - 0.5*(f(i,j-1)+f(i,j))*(jvel1!=0.0?1.0:0.0))/p->DYN[JP];
			  
    return grad;  
}

