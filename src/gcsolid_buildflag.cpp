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

void ghostcell::gcsolid_buildflag(lexer *p, fdm *a, int& cellcount)
{
    // Solid
    BASELOOP
    {
        if(a->solid(i,j,k)<0.0)
        p->flag4[IJK]=SOLID_FLAG;
			

        if(a->solid(i,j,k)>=0.0)
        p->flag4[IJK]=WATER_FLAG;
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
            if(p->flag4[Im1JK]<=SOLID_FLAG
            && p->flag4[Ip1JK]<=SOLID_FLAG)
            {
            p->flag4[IJK]=SOLID_FLAG;
            ++count;
            }
            
            if(p->j_dir==1)
            if(p->flag4[IJm1K]<=SOLID_FLAG
            && p->flag4[IJp1K]<=SOLID_FLAG)
            {
            p->flag4[IJK]=SOLID_FLAG;
            ++count;
            }

            if(p->k_dir==1)
            if(p->flag4[IJKm1]<=SOLID_FLAG
            && p->flag4[IJKp1]<=SOLID_FLAG)
            {
            p->flag4[IJK]=SOLID_FLAG;
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
    if(p->flag4[IJK]==SOLID_FLAG || (p->flag4[IJK]==WATER_FLAG && p->flag4[Ip1JK]==SOLID_FLAG))
	{
       if(p->flag4[IJK]==SOLID_FLAG) 
       p->flag1[IJK]=SOLID_FLAG;
       
       if(p->flag4[IJK]==WATER_FLAG && p->flag4[Ip1JK]==SOLID_FLAG)
       p->flag1[IJK]=SOLID_FLAG;
	}
	   
    if(p->flag4[IJK]==WATER_FLAG && p->flag4[Ip1JK]==WATER_FLAG)
    p->flag1[IJK]=WATER_FLAG;
    }
}

void ghostcell::gcsolid_velflag2(lexer *p, fdm *a, int& cellcount)
{
    count=0;
    BASELOOP
    {	
    if(p->flag4[IJK]==SOLID_FLAG || (p->flag4[IJK]==WATER_FLAG && p->flag4[IJp1K]==SOLID_FLAG))
	{
       if(p->flag4[IJK]==SOLID_FLAG) 
       p->flag2[IJK]=SOLID_FLAG;
       
       if(p->flag4[IJK]==WATER_FLAG && p->flag4[IJp1K]==SOLID_FLAG) 
       p->flag2[IJK]=SOLID_FLAG;
	}

    if(p->flag4[IJK]==WATER_FLAG && p->flag4[IJp1K]==WATER_FLAG)
    p->flag2[IJK]=WATER_FLAG;
    }

    cellcount=count;
}

void ghostcell::gcsolid_velflag3(lexer *p, fdm *a, int& cellcount)
{
    count=0;
    BASELOOP
    {
    if(p->flag4[IJK]==SOLID_FLAG || (p->flag4[IJK]==WATER_FLAG && p->flag4[IJKp1]==SOLID_FLAG))
	{
       if(p->flag4[IJK]==SOLID_FLAG) 
       p->flag3[IJK]=SOLID_FLAG;
       
       if(p->flag4[IJK]==WATER_FLAG && p->flag4[IJKp1]==SOLID_FLAG) 
       p->flag3[IJK]=SOLID_FLAG;
	}
	   
    if(p->flag4[IJK]==WATER_FLAG && p->flag4[IJKp1]==WATER_FLAG)
    p->flag3[IJK]=WATER_FLAG;
    }

    cellcount=count;
}
