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

void ghostcell::gcfb_buildflag(lexer *p, fdm *a, int **cellmem, int& cellcount)
{
    int cache;

    count=0;
    ALOOP // p-grid LOOP over i,j,k: flag4 = -17 and cellmem[3] = 2 if solid, flag4 = 10 and cellmem[3] = 1 if fluid
    {
    cache = p->flag4[IJK];

	

        if(a->fb(i,j,k)<0.0)
		{
        p->flag4[IJK]=FLT;
			
			if(cache>0)
			{
            cellmem[count][0]=i;
            cellmem[count][1]=j;
            cellmem[count][2]=k;
            cellmem[count][3]=2;
			
            ++count;
			}
		}

        if(a->fb(i,j,k)>=0.0)
        {
        p->flag4[IJK]=WATER;
		
			if(cache==FLT)
			{
            cellmem[count][0]=i;
            cellmem[count][1]=j;
            cellmem[count][2]=k;
            cellmem[count][3]=1;
            ++count;
			}
           
        }
		
		if(a->fb(i,j,k)>=0.0)
        p->flag4[IJK]=WATER;
    }
    cellcount=count;
    
    flagx(p,p->flag4);
    
	if(p->Y60==1)
    LOOP
    {
        if(p->flag4[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]==FLT
        && p->flag4[(i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]==FLT)
        p->flag4[IJK]=FLT;

        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]==FLT
        && p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin]==FLT)
        p->flag4[IJK]=FLT;

        if(p->flag4[IJK-1]==FLT
        && p->flag4[IJKp1]==FLT)
        p->flag4[IJK]=FLT;
    }
    
    cellcount=0;
    LOOP
    ++cellcount;
}


void ghostcell::gcfb_velflag1(lexer *p, fdm *a, int **cellmem, int& cellcount)
{
    int nn,cache;

    count=0;
    UBASELOOP // U-grid loop over i,j,k: flag1 = -17 and cellmem[3] = 2 if solid, flag1 = 10 and cellmem[3] = 1 if fluid
    {
	cache = p->flag1[IJK];
		
    if(p->flag4[IJK]<0 
	|| (p->flag4[IJK]>0 && p->flag4[(i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]<0))
	{
       p->flag1[IJK]=FLT;
	   
		if(cache>0)
		{
            cellmem[count][0]=i;
            cellmem[count][1]=j;
            cellmem[count][2]=k;
            cellmem[count][3]=2;
            ++count;
		}
	}
	   

    if(p->flag4[IJK]>0 && p->flag4[(i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]>0)
    {
		p->flag1[IJK]=WATER;
		
			if(cache==FLT)
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


void ghostcell::gcfb_velflag2(lexer *p, fdm *a, int **cellmem, int& cellcount)
{
	int nn,cache;

    count=0;
    VBASELOOP
    {
	cache = p->flag2[IJK];
		
    if(p->flag4[IJK]<0 
	|| (p->flag4[IJK]>0 && p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin]<0))
	{
       p->flag2[IJK]=FLT;
	   
		if(cache>0)
		{
            cellmem[count][0]=i;
            cellmem[count][1]=j;
            cellmem[count][2]=k;
            cellmem[count][3]=2;
            ++count;
		}
	}
	   

    if(p->flag4[IJK]>0 && p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin]>0)
    {
		p->flag2[IJK]=10;
		
			if(cache==FLT)
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

void ghostcell::gcfb_velflag3(lexer *p, fdm *a, int **cellmem, int& cellcount)
{
     int nn,cache;

    count=0;
    WBASELOOP
    {
	cache = p->flag3[IJK];
		
    if(p->flag4[IJK]<0 
	|| (p->flag4[IJK]>0 && p->flag4[IJKp1]<0))
	{
       p->flag3[IJK]=FLT;
	   
		if(cache>0)
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
		p->flag3[IJK]=10;
		
			if(cache==FLT)
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

