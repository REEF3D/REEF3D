/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"

void ghostcell::gcsolid_gcb_remove(lexer *p, fdm *a)
{   
    GCB4 
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];

    p->gcb4[n][3]=-fabs(p->gcb4[n][3]);
    }
}

void ghostcell::gcsolid_gcb_seed(lexer *p, fdm *a)
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
	
    // then check gcb4 around topo
	count=p->gcb_fix;
	LOOP
    {   
        // Solid
        if(p->flag4[Im1JK]==SOLID)
        ++count;
	
        if(p->flag4[IJp1K]==SOLID)
        ++count;

        if(p->flag4[IJm1K]==SOLID)
        ++count;

        if(p->flag4[Ip1JK]==SOLID)
        ++count;

        if(p->flag4[IJKm1]==SOLID)
        ++count;

        if(p->flag4[IJKp1]==SOLID)
        ++count;
    }
	
	p->Iresize(p->gcb4,p->gcb4_count, count, 6, 6); 
	p->Dresize(p->gcd4,p->gcb4_count, count); 
	
	count=p->gcb_fix;
	
	LOOP
    {
        // Solid
        if(p->flag4[Im1JK]==SOLID)
        {
        p->gcb4[count][0]=i;
        p->gcb4[count][1]=j;
        p->gcb4[count][2]=k;
        p->gcb4[count][3]=1;
        p->gcb4[count][4]=22;
        ++count;
        }

        if(p->flag4[IJp1K]==SOLID)
        {
        p->gcb4[count][0]=i;
        p->gcb4[count][1]=j;
        p->gcb4[count][2]=k;
        p->gcb4[count][3]=2;
        p->gcb4[count][4]=22;
        ++count;
        }

        if(p->flag4[IJm1K]==SOLID)
        {
        p->gcb4[count][0]=i;
        p->gcb4[count][1]=j;
        p->gcb4[count][2]=k;
        p->gcb4[count][3]=3;
        p->gcb4[count][4]=22;
        ++count;
        }

        if(p->flag4[Ip1JK]==SOLID)
        {
        p->gcb4[count][0]=i;
        p->gcb4[count][1]=j;
        p->gcb4[count][2]=k;
        p->gcb4[count][3]=4;
        p->gcb4[count][4]=22;
        ++count;
        }

        if(p->flag4[IJKm1]==SOLID)
        {
        p->gcb4[count][0]=i;
        p->gcb4[count][1]=j;
        p->gcb4[count][2]=k;
        p->gcb4[count][3]=5;
        p->gcb4[count][4]=22;
        ++count;
        }

        if(p->flag4[IJKp1]==SOLID)
        {
        p->gcb4[count][0]=i;
        p->gcb4[count][1]=j;
        p->gcb4[count][2]=k;
        p->gcb4[count][3]=6;
        p->gcb4[count][4]=22;
        ++count;
        }
    }
    p->gcb4_count=p->gcb_solid=p->gcb_topo=count;
	
	//cout<<p->mpirank<<" GCB4: "<<count<<endl;
}


void ghostcell::gcsolid_gcb_iocs(lexer *p, fdm *a)
{
    int count1,count2;

    count1=0;
    count2=0;
    GC4LOOP
    {
        if(fabs(p->gcb4[n][4])==1)
        ++count1;

        if(fabs(p->gcb4[n][4])==2)
        ++count2;
    }
	
	p->Iresize(p->gcin,p->gcin_count, count1, 6, 6);
	p->Iresize(p->gcout,p->gcout_count, count2, 6, 6); 	
	
	count1=0;
    count2=0;
    GC4LOOP
    {
        if(p->gcb4[n][4]==1)
        {
        p->gcin[count1][0]=p->gcb4[n][0];
        p->gcin[count1][1]=p->gcb4[n][1];
        p->gcin[count1][2]=p->gcb4[n][2];
        p->gcin[count1][3]=p->gcb4[n][3];
        p->gcin[count1][5]=p->gcb4[n][5];
        ++count1;
        }

        if(p->gcb4[n][4]==2)
        {
        p->gcout[count2][0]=p->gcb4[n][0];
        p->gcout[count2][1]=p->gcb4[n][1];
        p->gcout[count2][2]=p->gcb4[n][2];
        p->gcout[count2][3]=p->gcb4[n][3];
        p->gcout[count2][5]=p->gcb4[n][5];
        ++count2;
        }
    }
	
    p->gcin_count=count1;
    p->gcout_count=count2;
}



