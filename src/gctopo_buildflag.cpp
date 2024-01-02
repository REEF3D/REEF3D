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
#include"fdm.h"

void ghostcell::gcb_buildflag(lexer *p, fdm *a, int **cellmem, int& cellcount)
{
	int cache;
    
    //
    BASELOOP
    {
        if(a->solid(i,j,k)<0.0 || a->topo(i,j,k)<0.0)
        p->flag[IJK]=-1;
			

        if(a->solid(i,j,k)>=0.0 && a->topo(i,j,k)>=0.0)
        p->flag[IJK]=1;
    }
    
    // Topo
    count=0;
    SOLIDLOOP
    {
    cache = p->flag4[IJK];

        if(a->topo(i,j,k)<0.0)
		{
        p->flag4[IJK]=TOPO;
			
			if(cache>0)
			{
            cellmem[count][0]=i;
            cellmem[count][1]=j;
            cellmem[count][2]=k;
            cellmem[count][3]=2;
            ++count;
			}
		}

        if(a->topo(i,j,k)>=0.0)
        {
        p->flag4[IJK]=WATER;
		
			if(cache==TOPO)
			{
            cellmem[count][0]=i;
            cellmem[count][1]=j;
            cellmem[count][2]=k;
            cellmem[count][3]=1;
            ++count;
			}
        }
    }

    cellcount=count;
    
    flagx(p,p->flag);
    flagx(p,p->flag4);
    
	if(p->Y60==1)
    for(int qn=0; qn<100;++qn)
    {
        count=0;
        
        // check topo
        LOOP
        {
            if(p->i_dir==1)
            if(p->flag4[Im1JK]==TOPO
            && p->flag4[Ip1JK]==TOPO)
            {
            p->flag4[IJK]=TOPO;
            ++count;
            }
            
            if(p->j_dir==1)
            if(p->flag4[IJm1K]==TOPO
            && p->flag4[IJp1K]==TOPO)
            {
            p->flag4[IJK]=TOPO;
            ++count;
            }

            if(p->k_dir==1)
            if(p->flag4[IJKm1]==TOPO
            && p->flag4[IJKp1]==TOPO)
            {
            p->flag4[IJK]=TOPO;
            ++count;
            }
        }
        
        // check topo/solid/object combinations
        LOOP
        {
            if(p->i_dir==1)
            if(p->flag4[Im1JK]<0
            && p->flag4[Ip1JK]<0)
            {
            p->flag4[IJK]=TOPO;
            ++count;
            }
            
            if(p->j_dir==1)
            if(p->flag4[IJm1K]<0
            && p->flag4[IJp1K]<0)
            {
            p->flag4[IJK]=TOPO;
            ++count;
            }

            if(p->k_dir==1)
            if(p->flag4[IJKm1]<0
            && p->flag4[IJKp1]<0)
            {
            p->flag4[IJK]=TOPO;
            ++count;
            }
        }
        
        count = globalisum(count);
        
        if(count==0)
        break;
    }
    
    cellcount=0;
    LOOP
    ++cellcount;
}

void ghostcell::gcb_velflag1(lexer *p, fdm *a, int **cellmem, int& cellcount)
{
    int nn,cache;
    
    count=0;
    SOLIDLOOP
    {
	cache = p->flag1[IJK];
		
    if(p->flag4[IJK]<0 
	||(p->flag4[IJK]>0 && p->flag4[Ip1JK]<0))
	{
       if(p->flag4[IJK]<0) 
       p->flag1[IJK]=p->flag4[IJK];
       
       if(p->flag4[IJK]>0 && p->flag4[Ip1JK]<0)
       p->flag1[IJK]=p->flag4[Ip1JK];
	   
       
		if(cache>0 && p->flag1[IJK]==TOPO)
		{
            cellmem[count][0]=i;
            cellmem[count][1]=j;
            cellmem[count][2]=k;
            cellmem[count][3]=2;
            ++count;
		}
	}
	   

    if(p->flag4[IJK]>0 && p->flag4[Ip1JK]>0)
    {
		p->flag1[IJK]=WATER;
		
			if(cache==TOPO)
			{
            cellmem[count][0]=i;
            cellmem[count][1]=j;
            cellmem[count][2]=k;
            cellmem[count][3]=1;
            ++count;
			}
    }

    }
    

    cellcount=count;
}

void ghostcell::gcb_velflag2(lexer *p, fdm *a, int **cellmem, int& cellcount)
{
   int nn,cache;

    count=0;
    SOLIDLOOP
    {
	cache = p->flag2[IJK];
		
    if(p->flag4[IJK]<0 
	|| (p->flag4[IJK]>0 && p->flag4[IJp1K]<0))
	{
       if(p->flag4[IJK]<0) 
       p->flag2[IJK]=p->flag4[IJK];
       
       if(p->flag4[IJK]>0 && p->flag4[IJp1K]<0) 
       p->flag2[IJK]=p->flag4[IJp1K];
	   
		if(cache>0 && p->flag2[IJK]==TOPO)
		{
            cellmem[count][0]=i;
            cellmem[count][1]=j;
            cellmem[count][2]=k;
            cellmem[count][3]=2;
            ++count;
		}
	}
	   

    if(p->flag4[IJK]>0 && p->flag4[IJp1K]>0)
    {
		p->flag2[IJK]=WATER;
		
			if(cache==TOPO)
			{
            cellmem[count][0]=i;
            cellmem[count][1]=j;
            cellmem[count][2]=k;
            cellmem[count][3]=1;
            ++count;
			}
    }

    }

    cellcount=count;
}

void ghostcell::gcb_velflag3(lexer *p, fdm *a, int **cellmem, int& cellcount)
{
     int nn,cache;

    count=0;
    SOLIDLOOP
    {
	cache = p->flag3[IJK];
		
    if(p->flag4[IJK]<0 
	|| (p->flag4[IJK]>0 && p->flag4[IJKp1]<0))
	{
       if(p->flag4[IJK]<0) 
       p->flag3[IJK]=p->flag4[IJK];
       
       if(p->flag4[IJK]>0 && p->flag4[IJKp1]<0) 
       p->flag3[IJK]=p->flag4[IJKp1];
	   
		if(cache>0 && p->flag3[IJK]==TOPO)
		{
            cellmem[count][0]=i;
            cellmem[count][1]=j;
            cellmem[count][2]=k;
            cellmem[count][3]=2;
            ++count;
		}
	}
	   
    if(p->flag4[IJK]>0 && p->flag4[IJKp1]>0)
    {
		p->flag3[IJK]=WATER;
		
			if(cache==TOPO)
			{
            cellmem[count][0]=i;
            cellmem[count][1]=j;
            cellmem[count][2]=k;
            cellmem[count][3]=1;
            ++count;
			}
    }

    }

    cellcount=count;
}

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
        p->flag1[Im2JK] =-3;
        }
        
        if(p->gcb1[n][4]==2)
        {
        i=p->gcb1[n][0];
        j=p->gcb1[n][1];
        k=p->gcb1[n][2];
        
        p->flag1[Ip1JK] =-4;
        p->flag1[Ip2JK] =-4;
        p->flag1[Ip2JK] =-4;
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
        p->flag2[Im2JK] =-3;
        }
        
        if(p->gcb2[n][4]==2)
        {
        i=p->gcb2[n][0];
        j=p->gcb2[n][1];
        k=p->gcb2[n][2];
        
        p->flag2[Ip1JK] =-4;
        p->flag2[Ip2JK] =-4;
        p->flag2[Ip2JK] =-4;
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
        p->flag3[Im2JK] =-3;
        }
        
        if(p->gcb3[n][4]==2)
        {
        i=p->gcb3[n][0];
        j=p->gcb3[n][1];
        k=p->gcb3[n][2];
        
        p->flag3[Ip1JK] =-4;
        p->flag3[Ip2JK] =-4;
        p->flag3[Ip2JK] =-4;
        }    
    }
}
