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

void ghostcell::gcfb_seed(lexer *p, fdm *a)
{
    // reinstate wall, if reactive
    GCB4 
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];

    if(p->flag4[IJK]>0)
    p->gcb4[n][3]=fabs(p->gcb4[n][3]);
    
    if(p->flag4[IJK]<0)
    p->gcb4[n][3]=-fabs(p->gcb4[n][3]);
    }
    
    // then check gcb4 around fb
	count=p->gcb_topo;
	LOOP
    {
        if(p->flag4[Im1JK]==FLT)
        ++count;
	
        if(p->flag4[IJm1K]==FLT)
        ++count;

        if(p->flag4[IJp1K]==FLT)
        ++count;

        if(p->flag4[Ip1JK]==FLT)
        ++count;

        if(p->flag4[IJKm1]==FLT)
        ++count;

        if(p->flag4[IJKp1]==FLT)
        ++count;
    }
	
    //if(p->gcb4_count!=count)
	//cout<<p->mpirank<<" old gcb4: "<<p->gcb4_count<<" new gcb4: "<<count<<endl;
	
	p->Iresize(p->gcb4,p->gcb4_count, count, 6, 6); 
	p->Dresize(p->gcd4,p->gcb4_count, count); 
	
	p->gcb4_count=p->gcb_fb=count;
	
	
	count=p->gcb_topo;
	LOOP
    {
        if(p->flag4[Im1JK]==FLT)
        {
        p->gcb4[count][0]=i;
        p->gcb4[count][1]=j;
        p->gcb4[count][2]=k;
        p->gcb4[count][3]=1;
        p->gcb4[count][4]=41;
        ++count;
        }

        if(p->flag4[IJm1K]==FLT)
        {
        p->gcb4[count][0]=i;
        p->gcb4[count][1]=j;
        p->gcb4[count][2]=k;
        p->gcb4[count][3]=3;
        p->gcb4[count][4]=41;
        ++count;
        }

        if(p->flag4[IJp1K]==FLT)
        {
        p->gcb4[count][0]=i;
        p->gcb4[count][1]=j;
        p->gcb4[count][2]=k;
        p->gcb4[count][3]=2;
        p->gcb4[count][4]=41;
        ++count;
        }

        if(p->flag4[Ip1JK]==FLT)
        {
        p->gcb4[count][0]=i;
        p->gcb4[count][1]=j;
        p->gcb4[count][2]=k;
        p->gcb4[count][3]=4;
        p->gcb4[count][4]=41;
        ++count;
        }

        if(p->flag4[IJKm1]==FLT)
        {
        p->gcb4[count][0]=i;
        p->gcb4[count][1]=j;
        p->gcb4[count][2]=k;
        p->gcb4[count][3]=5;
        p->gcb4[count][4]=41;
        ++count;
        }

        if(p->flag4[IJKp1]==FLT)
        {
        p->gcb4[count][0]=i;
        p->gcb4[count][1]=j;
        p->gcb4[count][2]=k;
        p->gcb4[count][3]=6;
        p->gcb4[count][4]=41;
        ++count;
        }
    }
    
    
    //if(p->mpirank==0)
	//cout<<p->mpirank<<" 2. old gcb4: "<<p->gcb4_count<<" new gcb4: "<<count<<endl;
}
    

