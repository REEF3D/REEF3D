/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

void part::xchange_fill(lexer *p, ghostcell *pgc, int mode, double *F)
{
    for(q=0;q<6;++q)
    {
    sendcount[q]=0;
    recvcount[q]=0;
    }
    
    // find particles for xchange
    for(n=0;n<index;++n)
    if(Flag[n]>0)
    {
        if(mode==1)
        {
        i=p->posc_i(X[n]);
        j=p->posc_j(Y[n]);
        k=p->posc_k(Z[n]);
        }
        
        if(mode==2)
        {
        i=p->posc_i(XRK1[n]);
        j=p->posc_j(YRK1[n]);
        k=p->posc_k(ZRK1[n]);
        }
    
    
    if(p->flag5[IJK]==-1)
    {
    send[0][sendcount[0]] = F[n];
    ++sendcount[0];
    }
    
    if(p->flag5[IJK]==-2)
    {
    send[1][sendcount[1]] = F[n];
    ++sendcount[1];
    }

    if(p->flag5[IJK]==-3)
    {
    send[2][sendcount[2]] = F[n];
    ++sendcount[2];
    }  
    
    if(p->flag5[IJK]==-4)
    {
    send[3][sendcount[3]] = F[n];
    ++sendcount[3];
    }
    
    if(p->flag5[IJK]==-5)
    {
    send[4][sendcount[4]] = F[n];
    ++sendcount[4];
    } 
    
    if(p->flag5[IJK]==-6)
    {
    send[5][sendcount[5]] = F[n];
    ++sendcount[5];
    }
    }

}

