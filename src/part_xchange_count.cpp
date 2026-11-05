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
#include<limits>

void part::xchange_count(lexer *p, ghostcell *pgc, int mode)
{
    index_empty0 = index_empty;

    for(q=0;q<6;++q)
    {
        sendnum[q]=0;
        recvnum[q]=0;
    }

    if(xdest.size()<size_t(index))
    xdest.resize(index);
    
    if(xkx.size()<size_t(index))
    {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        xkx.resize(index,nan);  xky.resize(index,nan);  xkz.resize(index,nan);
        xki.resize(index,0);    xkj.resize(index,0);    xkk.resize(index,0);
    }

    // count particles for xchange; the target side is kept for xchange_sendid, the cell indices of a
    // position that did not change since the last call are not searched again (the parcels at rest)
    for(n=0;n<index;++n)
    {
        xdest[n] = 0;
        
        if(Flag[n]!=ACTIVE)
        continue;
        
        double xs = mode==1 ? XRK1[n] : X[n];
        double ys = mode==1 ? YRK1[n] : Y[n];
        double zs = mode==1 ? ZRK1[n] : Z[n];
        
        if(xs==xkx[n] && ys==xky[n] && zs==xkz[n])
        {
            i=xki[n];
            j=xkj[n];
            k=xkk[n];
        }
        else
        {
            i=p->posc_i(xs);
            j=p->posc_j(ys);
            k=p->posc_k(zs);
            
            xkx[n]=xs; xky[n]=ys; xkz[n]=zs;
            xki[n]=i;  xkj[n]=j;  xkk[n]=k;
        }

        int f = p->flag5[IJK];
        
        if(f<=-1 && f>=-6)
        {
            xdest[n] = -f;
            ++sendnum[-f-1];
        }
    }
}
