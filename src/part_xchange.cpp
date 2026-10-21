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

void part::xchange(lexer *p, ghostcell *pgc, slice &bedch, int mode)
{
    // count send/recv
    xchange_count(p,pgc,mode);
    pgc->gcpartnum(sendnum,recvnum);
    
    // nothing leaves any subdomain (the usual case for a resting bed): no transfer
    int nsend=0;
    for(q=0;q<6;++q)
    nsend += sendnum[q];
    
    if(pgc->globalisum(nsend)==0)
    return;

    // check send/recv array size: the parcel arrays go packed in one exchange
    const int NA = 19;
    xchange_resize(p,pgc,NA);

    // sendid
    xchange_sendid(p,pgc,mode);

    // the parcel arrays (after the resize, which may move them)
    double *arr[NA] = {U,V,W, URK1,VRK1,WRK1, X,Y,Z, XRK1,YRK1,ZRK1, Uf,Vf,Wf, D,RO, Test,Hop};

    for(n=0;n<6;++n)
    for(q=0;q<sendnum[n];++q)
    {
        const int m = sendid[n][q];
        double *b = &send[n][size_t(q)*NA];

        for(int a=0;a<NA;++a)
        b[a] = arr[a][m];
    }

    int sn[6], rn[6];
    for(q=0;q<6;++q)
    {
        sn[q] = NA*sendnum[q];
        rn[q] = NA*recvnum[q];
    }

    pgc->gcpartx(sn,rn,send,recv);

    // received parcels into the free slots, Empty[0..index_empty-1] holds them
    index_empty = index_empty0;

    for(n=0;n<6;++n)
    for(q=0;q<recvnum[n];++q)
    {
        --index_empty;
        const int m = Empty[index_empty];
        const double *b = &recv[n][size_t(q)*NA];

        for(int a=0;a<NA;++a)
        arr[a][m] = b[a];
    }

    // Flag
    xchange_fillback_flag(p,pgc,bedch,mode);
    xchange_fill_flag(p,pgc,mode);
}
