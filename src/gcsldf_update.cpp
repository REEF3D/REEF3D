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
#include"fieldint4.h"

void ghostcell::gcsldf_update(lexer *p)
{
    count=0;
    
    k=p->knoz-1;
    SLICELOOP4
    if(p->DF[IJK]==1)
    {
        if(p->DF[Im1JK]<0)
        ++count;
        
        if(p->DF[Ip1JK]<0)
        ++count;
        
        if(p->DF[IJm1K]<0)
        ++count;
        
        if(p->DF[IJp1K]<0)
        ++count;
    }
    
    
    if(p->gcsldf4_count!=count)
    {
    p->Iresize(p->gcsldf4,p->gcsldf4_count,count,6,6);
    
    p->gcsldf4_count=count;
    }
    
    
    
    
    // assign gcsldf entries
    
    count=0;
    
    k=p->knoz-1;
    SLICELOOP4
    if(p->DF[IJK]==1)
    {
        if(p->DF[Im1JK]<0)
        {
        p->gcsldf4[count][0]=i;
        p->gcsldf4[count][1]=j;
        p->gcsldf4[count][3]=1;
        p->gcsldf4[count][4]=48;
        ++count;
        }
        
        if(p->DF[Ip1JK]<0)
        {
        p->gcsldf4[count][0]=i;
        p->gcsldf4[count][1]=j;
        p->gcsldf4[count][3]=4;      
        p->gcsldf4[count][4]=48;
        ++count;
        }
        
        if(p->DF[IJm1K]<0)
        {
        p->gcsldf4[count][0]=i;
        p->gcsldf4[count][1]=j;
        p->gcsldf4[count][3]=3;
        p->gcsldf4[count][4]=48;
        ++count;
        }
        
        if(p->DF[IJp1K]<0)
        {
        p->gcsldf4[count][0]=i;
        p->gcsldf4[count][1]=j;
        p->gcsldf4[count][3]=2;
        p->gcsldf4[count][4]=48;
        ++count;
        }
    }
    
    
}