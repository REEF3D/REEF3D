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

#include"part.h"
#include"lexer.h"
#include"ghostcell.h"

void part::xchange_resize(lexer *p, ghostcell *pgc, int na)
{
    // check send / recv size: sendid holds a parcel per entry, send and recv na values per parcel
    maxnum=0;

    for(q=0;q<6;++q)
    {
        maxnum = MAX(maxnum,sendnum[q]);
        maxnum = MAX(maxnum,recvnum[q]);
    }

    if(maxnum>capacity_para)
    {
        p->Iresize(sendid,6,6,capacity_para,maxnum);

        capacity_para = maxnum;
    }

    if(na*maxnum>capacity_buf)
    {
        p->Dresize(send,6,6,capacity_buf,na*maxnum);
        p->Dresize(recv,6,6,capacity_buf,na*maxnum);

        capacity_buf = na*maxnum;
    }

    // check arrays
    maxnum=0;

    for(q=0;q<6;++q)
        maxnum += recvnum[q];

    int diff = index_empty - maxnum;

    if(diff<0)
    {
        resize(p,capacity + 2*fabs(diff));
        index_empty0 = index_empty;
    }
}
