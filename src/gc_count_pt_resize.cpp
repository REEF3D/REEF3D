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

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"

void ghostcell::column_pt_resize(lexer* p, fdm* a)
{
    int size=0;
    int safety=1000;
    
    
    size = column_pt4_count(p,a);
    
    if(size>p->C4_size)
    {
    cout<<p->mpirank<<" CPT4 Resize: "<<p->C4_size<<" "<<size<<endl;
    size+=safety;
    a->C4.resize(p,p->C4_size,size);
    
    p->C4_size=size;
    }
    
    
    size = column_pt4a_count(p,a);
    
    if(size>p->C4a_size)
    {
    cout<<p->mpirank<<" CPT4a Resize: "<<p->C4a_size<<" "<<size<<endl;
    size+=safety;
    a->C4a.resize(p,p->C4a_size,size);
    
    p->C4a_size=size;
    }
    
    
    size = column_pt6_count(p,a);
    
    if(size>p->C6_size)
    {
    cout<<p->mpirank<<" CPT6 Resize: "<<p->C6_size<<" "<<size<<endl;
    size+=safety;
    a->C6.resize(p,p->C6_size,size);
    
    p->C6_size=size;
    }
   
}
