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

void part::xchange_sendid(lexer *p, ghostcell *pgc, int mode)
{
    for(q=0;q<6;++q)
    sendnum[q]=0;

    
    // count particles for xchange
    for(n=0;n<index;++n)
    if(Flag[n]==ACTIVE)
    {
            if(mode==1)
            {
            i=p->posc_i(XRK1[n]);
            j=p->posc_j(YRK1[n]);
            k=p->posc_k(ZRK1[n]);
            }
            
            if(mode==2)
            {
            i=p->posc_i(X[n]);
            j=p->posc_j(Y[n]);
            k=p->posc_k(Z[n]);
            }
        
        if(p->flag5[IJK]==-1)
        {
        sendid[0][sendnum[0]]=n;
        ++sendnum[0]; 
        }
        
        if(p->flag5[IJK]==-2)
        {
        sendid[1][sendnum[1]]=n;
        ++sendnum[1];
        }

        if(p->flag5[IJK]==-3)
        {
        sendid[2][sendnum[2]]=n;
        ++sendnum[2];  
        }
        
        if(p->flag5[IJK]==-4)
        {
        sendid[3][sendnum[3]]=n;
        ++sendnum[3]; 
        }
        
        if(p->flag5[IJK]==-5)
        {
        sendid[4][sendnum[4]]=n;
        ++sendnum[4]; 
        }
        
        if(p->flag5[IJK]==-6)
        {
        sendid[5][sendnum[5]]=n;
        ++sendnum[5]; 
        }
    }
}

