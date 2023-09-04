/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

void ghostcell::gcsolid_buildflag(lexer *p, fdm *a, int& cellcount)
{
    // Solid
    BASELOOP
    {
        if(a->solid(i,j,k)<0.0)
        p->flag4[IJK]=SOLID;
			

        if(a->solid(i,j,k)>=0.0)
        p->flag4[IJK]=WATER;
    }
    
    //
    BASELOOP
    {
        if(a->solid(i,j,k)<0.0 || a->topo(i,j,k)<0.0)
        p->flag[IJK]=-1;
			

        if(a->solid(i,j,k)>=0.0 && a->topo(i,j,k)>=0.0)
        p->flag[IJK]=1;
    }
    
    flagx(p,p->flag);
    flagx(p,p->flag4);

	if(p->Y60==1)
    for(int qn=0; qn<100;++qn)
    {
        count=0;
        // check solid
        LOOP
        {
            if(p->i_dir==1)
            if(p->flag4[Im1JK]<=SOLID
            && p->flag4[Ip1JK]<=SOLID)
            {
            p->flag4[IJK]=SOLID;
            ++count;
            }
            
            if(p->j_dir==1)
            if(p->flag4[IJm1K]<=SOLID
            && p->flag4[IJp1K]<=SOLID)
            {
            p->flag4[IJK]=SOLID;
            ++count;
            }

            if(p->k_dir==1)
            if(p->flag4[IJKm1]<=SOLID
            && p->flag4[IJKp1]<=SOLID)
            {
            p->flag4[IJK]=SOLID;
            ++count;
            }
            
            
        }
        
        count = globalisum(count);
            
            //if(p->mpirank==0)
            //cout<<p->mpirank<<"  Y60 count: "<<count<<endl;
            
    if(count==0)
    break;
    }
    
    
    cellcount=0;
    LOOP
    ++cellcount;
}

void ghostcell::gcsolid_velflag1(lexer *p, fdm *a, int& cellcount)
{
    count=0;
    BASELOOP
    {
    if(p->flag4[IJK]==SOLID || (p->flag4[IJK]==WATER && p->flag4[Ip1JK]==SOLID))
	{
       if(p->flag4[IJK]==SOLID) 
       p->flag1[IJK]=SOLID;
       
       if(p->flag4[IJK]==WATER && p->flag4[Ip1JK]==SOLID)
       p->flag1[IJK]=SOLID;
	}
	   
    if(p->flag4[IJK]==WATER && p->flag4[Ip1JK]==WATER)
    p->flag1[IJK]=WATER;
    }
}

void ghostcell::gcsolid_velflag2(lexer *p, fdm *a, int& cellcount)
{
    count=0;
    BASELOOP
    {	
    if(p->flag4[IJK]==SOLID || (p->flag4[IJK]==WATER && p->flag4[IJp1K]==SOLID))
	{
       if(p->flag4[IJK]==SOLID) 
       p->flag2[IJK]=SOLID;
       
       if(p->flag4[IJK]==WATER && p->flag4[IJp1K]==SOLID) 
       p->flag2[IJK]=SOLID;
	}

    if(p->flag4[IJK]==WATER && p->flag4[IJp1K]==WATER)
    p->flag2[IJK]=WATER;
    }

    cellcount=count;
}

void ghostcell::gcsolid_velflag3(lexer *p, fdm *a, int& cellcount)
{
    count=0;
    BASELOOP
    {
    if(p->flag4[IJK]==SOLID || (p->flag4[IJK]==WATER && p->flag4[IJKp1]==SOLID))
	{
       if(p->flag4[IJK]==SOLID) 
       p->flag3[IJK]=SOLID;
       
       if(p->flag4[IJK]==WATER && p->flag4[IJKp1]==SOLID) 
       p->flag3[IJK]=SOLID;
	}
	   
    if(p->flag4[IJK]==WATER && p->flag4[IJKp1]==WATER)
    p->flag3[IJK]=WATER;
    }

    cellcount=count;
}
