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

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"

void ghostcell::gcb_velflagio(lexer *p, fdm *a)
{
    GC1LOOP
    {
        if(p->gcb1[n][4]==1)
        {
        i=p->gcb1[n][0];
        j=p->gcb1[n][1];
        k=p->gcb1[n][2];
        
        p->flag1[Im1JK] =-3;
        p->flag1[Im2JK] =-3;
        p->flag1[Im3JK] =-3;
        }
        
        if(p->gcb1[n][4]==2)
        {
        i=p->gcb1[n][0];
        j=p->gcb1[n][1];
        k=p->gcb1[n][2];
        
        p->flag1[Ip1JK] =-4;
        p->flag1[Ip2JK] =-4;
        p->flag1[Ip3JK] =-4;
        }    
    }
    
    
    GC2LOOP
    {
        if(p->gcb2[n][4]==1)
        {
        i=p->gcb2[n][0];
        j=p->gcb2[n][1];
        k=p->gcb2[n][2];
        
        p->flag2[Im1JK] =-3;
        p->flag2[Im2JK] =-3;
        p->flag2[Im3JK] =-3;
        }
        
        if(p->gcb2[n][4]==2)
        {
        i=p->gcb2[n][0];
        j=p->gcb2[n][1];
        k=p->gcb2[n][2];
        
        p->flag2[Ip1JK] =-4;
        p->flag2[Ip2JK] =-4;
        p->flag2[Ip3JK] =-4;
        }    
    }
    
    
    
    GC3LOOP
    {
        if(p->gcb3[n][4]==1)
        {
        i=p->gcb3[n][0];
        j=p->gcb3[n][1];
        k=p->gcb3[n][2];
        
        p->flag3[Im1JK] =-3;
        p->flag3[Im2JK] =-3;
        p->flag3[Im3JK] =-3;
        }
        
        if(p->gcb3[n][4]==2)
        {
        i=p->gcb3[n][0];
        j=p->gcb3[n][1];
        k=p->gcb3[n][2];
        
        p->flag3[Ip1JK] =-4;
        p->flag3[Ip2JK] =-4;
        p->flag3[Ip3JK] =-4;
        }    
    }
}
