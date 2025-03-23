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

void ghostcell::gcdf_update(lexer *p, fdm *a)
{
    // FLAGSF
    LOOP
    {
    if((a->fb(i,j,k)>0.0 || p->X10==0) && (a->solid(i,j,k)>0.0 || p->solidread==0) && (a->topo(i,j,k)>0.0 || p->toporead==0))
    p->flagsf4[IJK]=1;
    
    if((a->fb(i,j,k)<0.0 && p->X10==1) || (a->solid(i,j,k)<0.0 && p->solidread==1) || (a->topo(i,j,k)<0.0 && p->toporead==1))
    p->flagsf4[IJK]=-1;
    }
    
    flagx(p,p->flagsf4);
    
    LOOP
    {
    p->flagsf1[IJK]=p->flagsf4[IJK];
    p->flagsf2[IJK]=p->flagsf4[IJK];
    p->flagsf3[IJK]=p->flagsf4[IJK];
    
    if(p->flagsf4[IJK]>0 && p->flagsf4[Ip1JK]<0)
    p->flagsf1[IJK]=-1;
    
    if(p->flagsf4[IJK]>0 && p->flagsf4[IJp1K]<0)
    p->flagsf2[IJK]=-1;
    
    if(p->flagsf4[IJK]>0 && p->flagsf4[IJKp1]<0)
    p->flagsf3[IJK]=-1;
    }
    
    flagx(p,p->flagsf1);
    flagx(p,p->flagsf2);
    flagx(p,p->flagsf3);
        
    // count gcdf entries
    count=0;
    
    
    // flag 
    LOOP
    if(p->flagsf4[IJK]>0)
    {
     
        if(p->flagsf4[Im1JK]<0)
        ++count;
        
        if(p->flagsf4[Ip1JK]<0)
        ++count;
        
        if(p->flagsf4[IJm1K]<0)
        ++count;
        
        if(p->flagsf4[IJp1K]<0)
        ++count;

        if(p->flagsf4[IJKm1]<0)
        ++count;
        
        if(p->flagsf4[IJKp1]<0)
        ++count;        
    }
    
    if(p->gcdf4_count!=count)
    {
    p->Iresize(p->gcdf4,p->gcdf4_count,count,6,6);
    
    p->gcdf4_count=count;
    }
    
    //cout<<p->mpirank<<" p->gcdf4_count: "<<p->gcdf4_count<<endl;
    
    // assign gcdf entries
    count=0;
    
    LOOP
    if(p->flagsf4[IJK]>0)
    {
        if(p->flagsf4[Im1JK]<0)
        {
        p->gcdf4[count][0]=i;
        p->gcdf4[count][1]=j;
        p->gcdf4[count][2]=k;
        p->gcdf4[count][3]=1;
        p->gcdf4[count][4]=48;
        ++count;
        }
        
        if(p->flagsf4[Ip1JK]<0)
        {
        p->gcdf4[count][0]=i;
        p->gcdf4[count][1]=j;
        p->gcdf4[count][2]=k;
        p->gcdf4[count][3]=4;
        p->gcdf4[count][4]=48;
        ++count;
        }
        
        if(p->flagsf4[IJm1K]<0)
        {
        p->gcdf4[count][0]=i;
        p->gcdf4[count][1]=j;
        p->gcdf4[count][2]=k;
        p->gcdf4[count][3]=3;
        p->gcdf4[count][4]=48;
        ++count;
        }
        
        if(p->flagsf4[IJp1K]<0)
        {
        p->gcdf4[count][0]=i;
        p->gcdf4[count][1]=j;
        p->gcdf4[count][2]=k;
        p->gcdf4[count][3]=2;
        p->gcdf4[count][4]=48;
        ++count;
        }

        if(p->flagsf4[IJKm1]<0)
        {
        p->gcdf4[count][0]=i;
        p->gcdf4[count][1]=j;
        p->gcdf4[count][2]=k;
        p->gcdf4[count][3]=5;
        p->gcdf4[count][4]=48;
        ++count;
        }
        
        if(p->flagsf4[IJKp1]<0)
        {
        p->gcdf4[count][0]=i;
        p->gcdf4[count][1]=j;
        p->gcdf4[count][2]=k;
        p->gcdf4[count][3]=6;
        p->gcdf4[count][4]=48;
        ++count;
        }       
    }
    
    fieldint4 cval(p);
    
    count=0;

    FLUIDLOOP
	{
    cval(i,j,k)=count;
    
    ++count;
	}
    

    GCDF4LOOP
    {
    i=p->gcdf4[n][0];
    j=p->gcdf4[n][1];
    k=p->gcdf4[n][2];
	p->gcdf4[n][5]=cval(i,j,k);
	}
    
    
    // -----------------------
    // flagsf1
    
    LOOP
    if(p->flagsf1[IJK]>0)
    {
     
        if(p->flagsf1[Im1JK]<0)
        ++count;
        
        if(p->flagsf1[Ip1JK]<0)
        ++count;
        
        if(p->flagsf1[IJm1K]<0)
        ++count;
        
        if(p->flagsf1[IJp1K]<0)
        ++count;

        if(p->flagsf1[IJKm1]<0)
        ++count;
        
        if(p->flagsf1[IJKp1]<0)
        ++count;        
    }
    
    if(p->gcdf1_count!=count)
    {
    p->Iresize(p->gcdf1,p->gcdf1_count,count,6,6);
    
    p->gcdf1_count=count;
    }
    
    //cout<<p->mpirank<<" p->gcdf1_count: "<<p->gcdf1_count<<endl;
    
    // assign gcdf entries
    count=0;
    
    LOOP
    if(p->flagsf1[IJK]>0)
    {
        if(p->flagsf1[Im1JK]<0)
        {
        p->gcdf1[count][0]=i;
        p->gcdf1[count][1]=j;
        p->gcdf1[count][2]=k;
        p->gcdf1[count][3]=1;
        p->gcdf1[count][4]=48;
        ++count;
        }
        
        if(p->flagsf1[Ip1JK]<0)
        {
        p->gcdf1[count][0]=i;
        p->gcdf1[count][1]=j;
        p->gcdf1[count][2]=k;
        p->gcdf1[count][3]=4;
        p->gcdf1[count][4]=48;
        ++count;
        }
        
        if(p->flagsf1[IJm1K]<0)
        {
        p->gcdf1[count][0]=i;
        p->gcdf1[count][1]=j;
        p->gcdf1[count][2]=k;
        p->gcdf1[count][3]=3;
        p->gcdf1[count][4]=48;
        ++count;
        }
        
        if(p->flagsf1[IJp1K]<0)
        {
        p->gcdf1[count][0]=i;
        p->gcdf1[count][1]=j;
        p->gcdf1[count][2]=k;
        p->gcdf1[count][3]=2;
        p->gcdf1[count][4]=48;
        ++count;
        }

        if(p->flagsf1[IJKm1]<0)
        {
        p->gcdf1[count][0]=i;
        p->gcdf1[count][1]=j;
        p->gcdf1[count][2]=k;
        p->gcdf1[count][3]=5;
        p->gcdf1[count][4]=48;
        ++count;
        }
        
        if(p->flagsf1[IJKp1]<0)
        {
        p->gcdf1[count][0]=i;
        p->gcdf1[count][1]=j;
        p->gcdf1[count][2]=k;
        p->gcdf1[count][3]=6;
        p->gcdf1[count][4]=48;
        ++count;
        }       
    }
    
    // -----------------------
    // flagsf2
     LOOP
    if(p->flagsf2[IJK]>0)
    {
     
        if(p->flagsf2[Im1JK]<0)
        ++count;
        
        if(p->flagsf2[Ip1JK]<0)
        ++count;
        
        if(p->flagsf2[IJm1K]<0)
        ++count;
        
        if(p->flagsf2[IJp1K]<0)
        ++count;

        if(p->flagsf2[IJKm1]<0)
        ++count;
        
        if(p->flagsf2[IJKp1]<0)
        ++count;        
    }
    
    if(p->gcdf2_count!=count)
    {
    p->Iresize(p->gcdf2,p->gcdf2_count,count,6,6);
    
    p->gcdf2_count=count;
    }
    
    //cout<<p->mpirank<<" p->gcdf2_count: "<<p->gcdf2_count<<endl;
    
    // assign gcdf entries
    count=0;
    
    LOOP
    if(p->flagsf2[IJK]>0)
    {
        if(p->flagsf2[Im1JK]<0)
        {
        p->gcdf2[count][0]=i;
        p->gcdf2[count][1]=j;
        p->gcdf2[count][2]=k;
        p->gcdf2[count][3]=1;
        p->gcdf2[count][4]=48;
        ++count;
        }
        
        if(p->flagsf2[Ip1JK]<0)
        {
        p->gcdf2[count][0]=i;
        p->gcdf2[count][1]=j;
        p->gcdf2[count][2]=k;
        p->gcdf2[count][3]=4;
        p->gcdf2[count][4]=48;
        ++count;
        }
        
        if(p->flagsf2[IJm1K]<0)
        {
        p->gcdf2[count][0]=i;
        p->gcdf2[count][1]=j;
        p->gcdf2[count][2]=k;
        p->gcdf2[count][3]=3;
        p->gcdf2[count][4]=48;
        ++count;
        }
        
        if(p->flagsf2[IJp1K]<0)
        {
        p->gcdf2[count][0]=i;
        p->gcdf2[count][1]=j;
        p->gcdf2[count][2]=k;
        p->gcdf2[count][3]=2;
        p->gcdf2[count][4]=48;
        ++count;
        }

        if(p->flagsf2[IJKm1]<0)
        {
        p->gcdf2[count][0]=i;
        p->gcdf2[count][1]=j;
        p->gcdf2[count][2]=k;
        p->gcdf2[count][3]=5;
        p->gcdf2[count][4]=48;
        ++count;
        }
        
        if(p->flagsf2[IJKp1]<0)
        {
        p->gcdf2[count][0]=i;
        p->gcdf2[count][1]=j;
        p->gcdf2[count][2]=k;
        p->gcdf2[count][3]=6;
        p->gcdf2[count][4]=48;
        ++count;
        }       
    }
    
    // -----------------------
    // flagsf3
    
    LOOP
    if(p->flagsf3[IJK]>0)
    {
     
        if(p->flagsf3[Im1JK]<0)
        ++count;
        
        if(p->flagsf3[Ip1JK]<0)
        ++count;
        
        if(p->flagsf3[IJm1K]<0)
        ++count;
        
        if(p->flagsf3[IJp1K]<0)
        ++count;

        if(p->flagsf3[IJKm1]<0)
        ++count;
        
        if(p->flagsf3[IJKp1]<0)
        ++count;        
    }
    
    if(p->gcdf3_count!=count)
    {
    p->Iresize(p->gcdf3,p->gcdf3_count,count,6,6);
    
    p->gcdf3_count=count;
    }
    
    //cout<<p->mpirank<<" p->gcdf3_count: "<<p->gcdf3_count<<endl;
    
    // assign gcdf entries
    count=0;
    
    LOOP
    if(p->flagsf3[IJK]>0)
    {
        if(p->flagsf3[Im1JK]<0)
        {
        p->gcdf3[count][0]=i;
        p->gcdf3[count][1]=j;
        p->gcdf3[count][2]=k;
        p->gcdf3[count][3]=1;
        p->gcdf3[count][4]=48;
        ++count;
        }
        
        if(p->flagsf3[Ip1JK]<0)
        {
        p->gcdf3[count][0]=i;
        p->gcdf3[count][1]=j;
        p->gcdf3[count][2]=k;
        p->gcdf3[count][3]=4;
        p->gcdf3[count][4]=48;
        ++count;
        }
        
        if(p->flagsf3[IJm1K]<0)
        {
        p->gcdf3[count][0]=i;
        p->gcdf3[count][1]=j;
        p->gcdf3[count][2]=k;
        p->gcdf3[count][3]=3;
        p->gcdf3[count][4]=48;
        ++count;
        }
        
        if(p->flagsf3[IJp1K]<0)
        {
        p->gcdf3[count][0]=i;
        p->gcdf3[count][1]=j;
        p->gcdf3[count][2]=k;
        p->gcdf3[count][3]=2;
        p->gcdf3[count][4]=48;
        ++count;
        }

        if(p->flagsf3[IJKm1]<0)
        {
        p->gcdf3[count][0]=i;
        p->gcdf3[count][1]=j;
        p->gcdf3[count][2]=k;
        p->gcdf3[count][3]=5;
        p->gcdf3[count][4]=48;
        ++count;
        }
        
        if(p->flagsf3[IJKp1]<0)
        {
        p->gcdf3[count][0]=i;
        p->gcdf3[count][1]=j;
        p->gcdf3[count][2]=k;
        p->gcdf3[count][3]=6;
        p->gcdf3[count][4]=48;
        ++count;
        }       
    }


}