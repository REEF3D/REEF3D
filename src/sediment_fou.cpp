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

#include"sediment_fou.h"
#include"lexer.h"
#include"slice.h"
#include"vec.h"
#include"fnpf_discrete_weights.h"

sediment_fou::sediment_fou(lexer* p) 
{
}

sediment_fou::~sediment_fou()
{
}

double sediment_fou::sx(lexer *p, slice &f, double ivel1, double ivel2)
{
    if(ivel1>=0.0)
    fu1 = f(i-1,j);
    
    if(ivel1<0.0)
    fu1 = f(i,j);
    
    if(ivel2>=0.0)
    fu2 = f(i,j);
    
    if(ivel2<0.0)
    fu2 = f(i+1,j);
        
        if(p->S31==1)
        grad = ((fu2*ivel2-fu1*ivel1)/p->DXN[IP]);
        
        // S 31 2: f is the flux component of the cells, qb*s_x. Face flux by flux splitting,
        // max(f_left,0) + min(f_right,0): every cell sends its own flux through the face it points to. The face
        // direction only closes the face (0 next to structures and walls). Taking the upwind cell by the sign of
        // the face velocity gave the whole flux of one cell or of the other where that velocity is zero up to
        // round-off (symmetry lines, stagnation points): symmetric cases lost their symmetry and the bed
        // depended on the number of ranks.
        if(p->S31==2)
        {
        fu1 = (ivel1!=0.0) ? MAX(f(i-1,j),0.0) + MIN(f(i,j),0.0) : 0.0;
        fu2 = (ivel2!=0.0) ? MAX(f(i,j),0.0) + MIN(f(i+1,j),0.0) : 0.0;
        
        grad = (fu2 - fu1)/p->DXN[IP];
        }
        
        // S 31 3: upwind value of qb, closed faces (zero face direction) carry none
        if(p->S31==3)
        grad = ((fu2*(ivel2!=0.0?1.0:0.0) - fu1*(ivel1!=0.0?1.0:0.0))/p->DXN[IP]);
        
    return grad;
}

double sediment_fou::sy(lexer *p, slice &f, double jvel1, double jvel2)
{
    if(jvel1>=0.0)
    fv1 = f(i,j-1);
    
    if(jvel1<0.0)
    fv1 = f(i,j);
    
    if(jvel2>=0.0)
    fv2 = f(i,j);
    
    if(jvel2<0.0)
    fv2 = f(i,j+1);
    
        if(p->S31==1)
        grad = ((fv2*jvel2-fv1*jvel1)/p->DYN[JP]);
        
        if(p->S31==2)
        {
        fv1 = (jvel1!=0.0) ? MAX(f(i,j-1),0.0) + MIN(f(i,j),0.0) : 0.0;
        fv2 = (jvel2!=0.0) ? MAX(f(i,j),0.0) + MIN(f(i,j+1),0.0) : 0.0;
        
        grad = (fv2 - fv1)/p->DYN[JP];
        }
        
        if(p->S31==3)
        grad = ((fv2*(jvel2!=0.0?1.0:0.0) - fv1*(jvel1!=0.0?1.0:0.0))/p->DYN[JP]);
			  
    return grad;  
}

