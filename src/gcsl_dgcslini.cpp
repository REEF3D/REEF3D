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

#include"ghostcell.h"
#include"lexer.h"
#include"sliceint5.h"

void ghostcell::dgcslini1(lexer* p)
{
    count=0;
    SLICELOOP1
    {
        if(p->flagslice1[(i-p->imin-1)*p->jmax + (j-p->jmin)]<0
        && p->flagslice1[(i-p->imin)*p->jmax + (j-p->jmin-1)]<0
        && p->flagslice1[(i-p->imin-1)*p->jmax + (j-p->jmin-1)]<0)
        {
        ++count;
        }

        if(p->flagslice1[(i-p->imin+1)*p->jmax + (j-p->jmin)]<0
        && p->flagslice1[(i-p->imin)*p->jmax + (j-p->jmin)-1]<0
        && p->flagslice1[(i-p->imin+1)*p->jmax + (j-p->jmin)-1]<0)
        {
        ++count;
        }
        
        if(p->flagslice1[(i-p->imin+1)*p->jmax + (j-p->jmin)]<0
        && p->flagslice1[(i-p->imin)*p->jmax + (j-p->jmin)+1]<0
        && p->flagslice1[(i-p->imin+1)*p->jmax + (j-p->jmin)+1]<0)
        {
        ++count;
        }
        
        if(p->flagslice1[(i-p->imin-1)*p->jmax + (j-p->jmin)]<0
        && p->flagslice1[(i-p->imin)*p->jmax + (j-p->jmin)+1]<0
        && p->flagslice1[(i-p->imin-1)*p->jmax + (j-p->jmin)+1]<0)
        {
        ++count;
        }
    }
    
    p->Iarray(p->dgcsl1,count,3);
	p->dgcsl1_count=count;
    
        //------------
        
    count=0;   
    SLICELOOP1
    {
        if(p->flagslice1[(i-p->imin-1)*p->jmax + (j-p->jmin)]<0
        && p->flagslice1[(i-p->imin)*p->jmax + (j-p->jmin-1)]<0
        && p->flagslice1[(i-p->imin-1)*p->jmax + (j-p->jmin-1)]<0)
        {
        p->dgcsl1[count][0]=i;
        p->dgcsl1[count][1]=j;
        p->dgcsl1[count][2]=1;
        ++count;
        }

        if(p->flagslice1[(i-p->imin+1)*p->jmax + (j-p->jmin)]<0
        && p->flagslice1[(i-p->imin)*p->jmax + (j-p->jmin)-1]<0
        && p->flagslice1[(i-p->imin+1)*p->jmax + (j-p->jmin)-1]<0)
        {
        p->dgcsl1[count][0]=i;
        p->dgcsl1[count][1]=j;
        p->dgcsl1[count][2]=2;
        ++count;
        }
        
        if(p->flagslice1[(i-p->imin+1)*p->jmax + (j-p->jmin)]<0
        && p->flagslice1[(i-p->imin)*p->jmax + (j-p->jmin)+1]<0
        && p->flagslice1[(i-p->imin+1)*p->jmax + (j-p->jmin)+1]<0)
        {
        p->dgcsl1[count][0]=i;
        p->dgcsl1[count][1]=j;
        p->dgcsl1[count][2]=3;
        ++count;
        }
        
        if(p->flagslice1[(i-p->imin-1)*p->jmax + (j-p->jmin)]<0
        && p->flagslice1[(i-p->imin)*p->jmax + (j-p->jmin)+1]<0
        && p->flagslice1[(i-p->imin-1)*p->jmax + (j-p->jmin)+1]<0)
        {
        p->dgcsl1[count][0]=i;
        p->dgcsl1[count][1]=j;
        p->dgcsl1[count][2]=4;
        ++count;
        }
    }
    
}

void ghostcell::dgcslini2(lexer* p)
{
   count=0;
    SLICELOOP2
    {
        if(p->flagslice2[(i-p->imin-1)*p->jmax + (j-p->jmin)]<0
        && p->flagslice2[(i-p->imin)*p->jmax + (j-p->jmin-1)]<0
        && p->flagslice2[(i-p->imin-1)*p->jmax + (j-p->jmin-1)]<0)
        {
        ++count;
        }

        if(p->flagslice2[(i-p->imin+1)*p->jmax + (j-p->jmin)]<0
        && p->flagslice2[(i-p->imin)*p->jmax + (j-p->jmin)-1]<0
        && p->flagslice2[(i-p->imin+1)*p->jmax + (j-p->jmin)-1]<0)
        {
        ++count;
        }
        
        if(p->flagslice2[(i-p->imin+1)*p->jmax + (j-p->jmin)]<0
        && p->flagslice2[(i-p->imin)*p->jmax + (j-p->jmin)+1]<0
        && p->flagslice2[(i-p->imin+1)*p->jmax + (j-p->jmin)+1]<0)
        {
        ++count;
        }
        
        if(p->flagslice2[(i-p->imin-1)*p->jmax + (j-p->jmin)]<0
        && p->flagslice2[(i-p->imin)*p->jmax + (j-p->jmin)+1]<0
        && p->flagslice2[(i-p->imin-1)*p->jmax + (j-p->jmin)+1]<0)
        {
        ++count;
        }
    }
    
    p->Iarray(p->dgcsl2,count,3);
	p->dgcsl2_count=count;
    
        //------------
        
    count=0;   
    SLICELOOP2
    {
        if(p->flagslice2[(i-p->imin-1)*p->jmax + (j-p->jmin)]<0
        && p->flagslice2[(i-p->imin)*p->jmax + (j-p->jmin-1)]<0
        && p->flagslice2[(i-p->imin-1)*p->jmax + (j-p->jmin-1)]<0)
        {
        p->dgcsl2[count][0]=i;
        p->dgcsl2[count][1]=j;
        p->dgcsl2[count][2]=1;
        ++count;
        }

        if(p->flagslice2[(i-p->imin+1)*p->jmax + (j-p->jmin)]<0
        && p->flagslice2[(i-p->imin)*p->jmax + (j-p->jmin)-1]<0
        && p->flagslice2[(i-p->imin+1)*p->jmax + (j-p->jmin)-1]<0)
        {
        p->dgcsl2[count][0]=i;
        p->dgcsl2[count][1]=j;
        p->dgcsl2[count][2]=2;
        ++count;
        }
        
        if(p->flagslice2[(i-p->imin+1)*p->jmax + (j-p->jmin)]<0
        && p->flagslice2[(i-p->imin)*p->jmax + (j-p->jmin)+1]<0
        && p->flagslice2[(i-p->imin+1)*p->jmax + (j-p->jmin)+1]<0)
        {
        p->dgcsl2[count][0]=i;
        p->dgcsl2[count][1]=j;
        p->dgcsl2[count][2]=3;
        ++count;
        }
        
        if(p->flagslice2[(i-p->imin-1)*p->jmax + (j-p->jmin)]<0
        && p->flagslice2[(i-p->imin)*p->jmax + (j-p->jmin)+1]<0
        && p->flagslice2[(i-p->imin-1)*p->jmax + (j-p->jmin)+1]<0)
        {
        p->dgcsl2[count][0]=i;
        p->dgcsl2[count][1]=j;
        p->dgcsl2[count][2]=4;
        ++count;
        }
    }
    		
}

void ghostcell::dgcslini3(lexer* p)
{
		
}

void ghostcell::dgcslini4(lexer* p)
{
    count=0;
    SLICELOOP4
    {
        if(p->flagslice4[(i-p->imin-1)*p->jmax + (j-p->jmin)]<0
        && p->flagslice4[(i-p->imin)*p->jmax + (j-p->jmin-1)]<0
        && p->flagslice4[(i-p->imin-1)*p->jmax + (j-p->jmin-1)]<0)
        {
        ++count;
        }

        if(p->flagslice4[(i-p->imin+1)*p->jmax + (j-p->jmin)]<0
        && p->flagslice4[(i-p->imin)*p->jmax + (j-p->jmin)-1]<0
        && p->flagslice4[(i-p->imin+1)*p->jmax + (j-p->jmin)-1]<0)
        {
        ++count;
        }
        
        if(p->flagslice4[(i-p->imin+1)*p->jmax + (j-p->jmin)]<0
        && p->flagslice4[(i-p->imin)*p->jmax + (j-p->jmin)+1]<0
        && p->flagslice4[(i-p->imin+1)*p->jmax + (j-p->jmin)+1]<0)
        {
        ++count;
        }
        
        if(p->flagslice4[(i-p->imin-1)*p->jmax + (j-p->jmin)]<0
        && p->flagslice4[(i-p->imin)*p->jmax + (j-p->jmin)+1]<0
        && p->flagslice4[(i-p->imin-1)*p->jmax + (j-p->jmin)+1]<0)
        {
        ++count;
        }
    }
    
    p->Iarray(p->dgcsl4,count,3);
	p->dgcsl4_count=count;
    
        //------------
        
    count=0;   
    SLICELOOP4
    {
        if(p->flagslice4[(i-p->imin-1)*p->jmax + (j-p->jmin)]<0
        && p->flagslice4[(i-p->imin)*p->jmax + (j-p->jmin-1)]<0
        && p->flagslice4[(i-p->imin-1)*p->jmax + (j-p->jmin-1)]<0)
        {
        p->dgcsl4[count][0]=i;
        p->dgcsl4[count][1]=j;
        p->dgcsl4[count][2]=1;
        ++count;
        }

        if(p->flagslice4[(i-p->imin+1)*p->jmax + (j-p->jmin)]<0
        && p->flagslice4[(i-p->imin)*p->jmax + (j-p->jmin)-1]<0
        && p->flagslice4[(i-p->imin+1)*p->jmax + (j-p->jmin)-1]<0)
        {
        p->dgcsl4[count][0]=i;
        p->dgcsl4[count][1]=j;
        p->dgcsl4[count][2]=2;
        ++count;
        }
        
        if(p->flagslice4[(i-p->imin+1)*p->jmax + (j-p->jmin)]<0
        && p->flagslice4[(i-p->imin)*p->jmax + (j-p->jmin)+1]<0
        && p->flagslice4[(i-p->imin+1)*p->jmax + (j-p->jmin)+1]<0)
        {
        p->dgcsl4[count][0]=i;
        p->dgcsl4[count][1]=j;
        p->dgcsl4[count][2]=3;
        ++count;
        }
        
        if(p->flagslice4[(i-p->imin-1)*p->jmax + (j-p->jmin)]<0
        && p->flagslice4[(i-p->imin)*p->jmax + (j-p->jmin)+1]<0
        && p->flagslice4[(i-p->imin-1)*p->jmax + (j-p->jmin)+1]<0)
        {
        p->dgcsl4[count][0]=i;
        p->dgcsl4[count][1]=j;
        p->dgcsl4[count][2]=4;
        ++count;
        }
    }
    
}
