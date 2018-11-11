/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include"lexer.h"
#include"fdm2D.h"

void ghostcell::column2D_pt_resize(lexer* p, fdm2D* b)
{
    int size=0;
    int safety=100;
    
    size = column2D_pt1_count(p,b);
    
    if(size>p->C1_2D_size)
    {
    cout<<p->mpirank<<" CPT2D_1 Resize: "<<p->C1_2D_size<<" "<<size<<endl;
    size+=safety;
    b->C1.resize(p,p->C4_size,size);
    
    p->C1_2D_size=size;
    }
    
    
    
    size = column2D_pt2_count(p,b);
    
    if(size>p->C2_2D_size)
    {
    cout<<p->mpirank<<" CPT2D_2 Resize: "<<p->C2_2D_size<<" "<<size<<endl;
    size+=safety;
    b->C2.resize(p,p->C2_2D_size,size);
    
    p->C2_2D_size=size;
    }
    
    
    
    size = column2D_pt4_count(p,b);
    
    if(size>p->C4_2D_size)
    {
    cout<<p->mpirank<<" CPT2D_4 Resize: "<<p->C4_2D_size<<" "<<size<<endl;
    size+=safety;
    b->C4.resize(p,p->C4_2D_size,size);
    
    p->C4_2D_size=size;
    }

}