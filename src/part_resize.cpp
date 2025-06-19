/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
Author: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"part.h"
#include"lexer.h"

void part::resize(lexer *p, int capacity_new)
{
    capacity_new+=1000;
    
    if(capacity_new>capacity)
    {
    p->Dresize(U,capacity,capacity_new);
    p->Dresize(V,capacity,capacity_new);
    p->Dresize(W,capacity,capacity_new);
    
    p->Dresize(URK1,capacity,capacity_new);
    p->Dresize(VRK1,capacity,capacity_new);
    p->Dresize(WRK1,capacity,capacity_new); 

    p->Dresize(X,capacity,capacity_new);
    p->Dresize(Y,capacity,capacity_new);
    p->Dresize(Z,capacity,capacity_new);
    
    p->Dresize(XRK1,capacity,capacity_new);
    p->Dresize(YRK1,capacity,capacity_new);
    p->Dresize(ZRK1,capacity,capacity_new);  
    
    p->Dresize(Uf,capacity,capacity_new);
    p->Dresize(Vf,capacity,capacity_new);
    p->Dresize(Wf,capacity,capacity_new);
    
    p->Dresize(D,capacity,capacity_new);
    p->Dresize(RO,capacity,capacity_new);
    
    p->Dresize(Test,capacity,capacity_new);
    
    p->Iresize(Flag,capacity,capacity_new);
    p->Iresize(Empty,capacity,capacity_new);
    
    // Flag update
    // add new cells to Empty list
    
    for(n=capacity;n<capacity_new;++n)
    {
    Flag[n] = -1;
    
    Empty[index_empty] = n;
    ++index_empty;
    }
    
    //cout<<p->mpirank<<" part.resize  CAPACITY_NEW: "<<capacity_new<<" CAPACITY: "<<capacity<<" index_empty: "<<index_empty<<endl;
    
    
    index=capacity=capacity_new;
    }
}