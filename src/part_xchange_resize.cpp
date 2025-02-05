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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"part.h"
#include"lexer.h"
#include"ghostcell.h"

void part::xchange_resize(lexer *p, ghostcell *pgc)
{
    // check send / recv size
    maxnum=0;
    
    for(q=0;q<6;++q)
    {
    maxnum = MAX(maxnum,sendnum[q]);
    maxnum = MAX(maxnum,recvnum[q]);
    }
    
    if(maxnum>capacity_para)
    {
    p->Dresize(send,6,6,capacity_para,maxnum);
    p->Dresize(recv,6,6,capacity_para,maxnum);
    
    p->Iresize(sendid,6,6,capacity_para,maxnum);
    
    capacity_para = maxnum;
    }
    
    // check arrays
    maxnum=0;
    
    for(q=0;q<6;++q)
    {
    maxnum += recvnum[q];
    }
    
    int diff = index_empty - maxnum;
    
    //cout<<p->mpirank<<" index_empty: "<<index_empty<<" maxnum: "<<maxnum<<" diff: "<<diff<<endl;

    if(diff<0)
    {
    resize(p,capacity + 2*fabs(diff));
    index_empty0 = index_empty;
    }
}


