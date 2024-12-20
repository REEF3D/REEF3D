/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is boundary of REEF3D.

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
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"boundary.h"
#include"lexer.h"

boundary::boundary(lexer *p, ghostcell *pgc)
{	
    capacity=1;
    index=1;
    index_empty=capacity;


    // 
    p->Iarray(iloc,capacity);
    p->Iarray(jloc,capacity);
    p->Iarray(kloc,capacity);
    
    p->Iarray(cellside,capacity);
    p->Iarray(bc_type,capacity);
    
    p->Darray(ks,capacity);
  
}

boundary::~boundary()
{
    delete[] iloc;
    iloc=nullptr;
    
    delete[] jloc;
    jloc=nullptr;
    
    delete[] kloc;
    kloc=nullptr;
    
    delete[] cellside;
    cellside=nullptr;
    
    delete[] bc_type;
    bc_type=nullptr;
    
    delete[] ks;
    ks=nullptr;
    
}


