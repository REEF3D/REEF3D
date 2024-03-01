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

#include"mgc3.h"
#include"lexer.h"

void mgc3::fillgcb(lexer *p)
{
    int q,n;
	
	p->Iarray(p->fgc,imax*jmax*kmax,6);

//  ------------
	
	if(p->gcb3_count!=p->gcb4_count)
	{
	p->Iresize(p->gcb3,p->gcb3_count, p->gcb4_count, 6, 6); 
	p->Dresize(p->gcd3,p->gcb3_count, p->gcb4_count);
	
	p->gcb3_count=p->gcb4_count;
	}
	
	QGCB4
	{
	for(n=0;n<5;++n)
	p->gcb3[q][n]=p->gcb4[q][n];

	if(p->gcb3[q][3]==5 || p->gcb3[q][3]==6)
	p->gcd3[q]=p->gcd4[q];

	if(p->gcb3[q][3]!=5 && p->gcb3[q][3]!=6)
	p->gcd3[q]=p->gcd4[q];
	}

    QGC3LOOP
	{
	    i=p->gcb3[q][0];
		j=p->gcb3[q][1];
		k=p->gcb3[q][2];
        
        if(p->gcb3[q][3]==5 || p->gcb3[q][3]==6)
        p->gcd3[q] += 0.5*p->DZP[KP];

		p->fgc[IJK][p->gcb3[q][3]-1]=1;
	}
	
	
	
	QGC3LOOP
	{
	    i=p->gcb3[q][0];
		j=p->gcb3[q][1];
		k=p->gcb3[q][2];

            if(p->gcb3[q][3]==6 && (p->periodic3!=1 || k+p->origin_k<p->gknoz-1))
            p->gcb3[q][2]-=1;
	}
    
    QGC3LOOP
	{
	    i=p->gcb3[q][0];
		j=p->gcb3[q][1];
		k=p->gcb3[q][2];

            if(p->gcb3[q][3]!=6 && p->fgc[IJK][5]==1 && (p->periodic3!=1 || k+p->origin_k<p->gknoz-1))
            p->gcb3[q][3]=-fabs(p->gcb3[q][3]);
	}
}

void mgc3::extragcb(lexer *p)
{
    int count;
    int q,n;

	for(n=0;n<imax*jmax*kmax;++n)
	for(q=0;q<6;++q)
	p->fgc[n][q]=0;

    QGC3LOOP
	{
	    i=p->gcb3[q][0];
		j=p->gcb3[q][1];
		k=p->gcb3[q][2];

		p->fgc[IJK][p->gcb3[q][3]-1]=1;
	}

    count=p->gcb3_count;
	
	// Find number of extragcb
	
	WLOOP
    {
        if(p->flag3[Im1JK]<0)
        if(p->fgc[IJK][0]==0)
        ++count;

        if(p->flag3[IJm1K]<0)
        if(p->fgc[IJK][2]==0)
        ++count;

        if(p->flag3[IJp1K]<0)
        if(p->fgc[IJK][1]==0)
        ++count;

        if(p->flag3[Ip1JK]<0)
        if(p->fgc[IJK][3]==0)
        ++count;

        if(p->flag3[IJKm1]<0)
        if(p->fgc[IJK][4]==0)
        ++count;

        if(p->flag3[IJKp1]<0)
        if(p->fgc[IJK][5]==0)
        ++count;
    }
    
    // extra parax + gcb
    for(q=0;q<p->gcpara1_count;++q)
    {
    i=p->gcpara1[q][0];
    j=p->gcpara1[q][1];
    k=p->gcpara1[q][2];
        
        if(p->flag3[IJK]==0 && p->flag3[Im1JK]==1)
        ++count;
    }

    for(q=0;q<p->gcpara3_count;++q)
    {
    i=p->gcpara3[q][0];
    j=p->gcpara3[q][1];
    k=p->gcpara3[q][2];
        
        if(p->flag3[IJK]==0 && p->flag3[IJm1K]==1)
        ++count;
    }

	for(q=0;q<p->gcpara5_count;++q)
	{
    i=p->gcpara5[q][0];
    j=p->gcpara5[q][1];
    k=p->gcpara5[q][2];
        
        if(p->flag3[IJK]==0 && p->flag3[IJKm1]==1)
        ++count;
	}

	for(q=0;q<p->gcpara4_count;++q)
	{
    i=p->gcpara4[q][0]+1;
    j=p->gcpara4[q][1];
    k=p->gcpara4[q][2];
        
        if(p->flag3[IJK]==0 && p->flag3[Ip1JK]==1)
        ++count;
	}

	for(q=0;q<p->gcpara2_count;++q)
	{
    i=p->gcpara2[q][0];
    j=p->gcpara2[q][1]+1;
    k=p->gcpara2[q][2];
        
        if(p->flag3[IJK]==0 && p->flag3[IJp1K]==1)
        ++count;
	}

	for(q=0;q<p->gcpara6_count;++q)
	{
	i=p->gcpara6[q][0];
    j=p->gcpara6[q][1];
    k=p->gcpara6[q][2]+1;
        
        if(p->flag3[IJK]==0 && p->flag3[IJKp1]==1)
        ++count;
	}
    
    // extra para + gcb 

	if(p->gcb3_count!=count)
	{
	p->Iresize(p->gcb3,p->gcb3_count, count, 6, 6); 
	p->Dresize(p->gcd3,p->gcb3_count, count); 
	}
	
// Store extragcb in gcb vec	
	count=p->gcb3_count;

    WLOOP
    {
        if(p->flag3[Im1JK]<0)
        if(p->fgc[IJK][0]==0)
        {
        p->gcb3[count][0]=i;
        p->gcb3[count][1]=j;
        p->gcb3[count][2]=k;
        p->gcb3[count][3]=1;

        if(p->flag3[Im1JK]==TOPO)
        p->gcb3[count][4]=5;

        if(p->flag3[Im1JK]==OBJ)
        p->gcb3[count][4]=21;
        
        if(p->flag3[Im1JK]==SOLID)
        p->gcb3[count][4]=22;
		
		if(p->flag3[Im1JK]==FLT)
        p->gcb3[count][4]=41;
        ++count;
        }

        if(p->flag3[IJm1K]<0)
        if(p->fgc[IJK][2]==0)
        {
        p->gcb3[count][0]=i;
        p->gcb3[count][1]=j;
        p->gcb3[count][2]=k;
        p->gcb3[count][3]=3;

        if(p->flag3[IJm1K]==TOPO)
        p->gcb3[count][4]=5;

        if(p->flag3[IJm1K]==OBJ)
        p->gcb3[count][4]=21;
        
        if(p->flag3[IJm1K]==SOLID)
        p->gcb3[count][4]=22;
		
		if(p->flag3[IJm1K]==FLT)
        p->gcb3[count][4]=41;
        ++count;
        }

        if(p->flag3[IJp1K]<0)
        if(p->fgc[IJK][1]==0)
        {
        p->gcb3[count][0]=i;
        p->gcb3[count][1]=j;
        p->gcb3[count][2]=k;
        p->gcb3[count][3]=2;

        if(p->flag3[IJp1K]==TOPO)
        p->gcb3[count][4]=5;

        if(p->flag3[IJp1K]==OBJ)
        p->gcb3[count][4]=21;
        
        if(p->flag3[IJm1K]==SOLID)
        p->gcb3[count][4]=22;
		
		if(p->flag3[IJp1K]==FLT)
        p->gcb3[count][4]=41;
        ++count;
        }

        if(p->flag3[Ip1JK]<0)
        if(p->fgc[IJK][3]==0)
        {
        p->gcb3[count][0]=i;
        p->gcb3[count][1]=j;
        p->gcb3[count][2]=k;
        p->gcb3[count][3]=4;

        if(p->flag3[Ip1JK]==TOPO)
        p->gcb3[count][4]=5;

        if(p->flag3[Ip1JK]==OBJ)
        p->gcb3[count][4]=21;
        
        if(p->flag3[Ip1JK]==SOLID)
        p->gcb3[count][4]=22;
		
		if(p->flag3[Ip1JK]==FLT)
        p->gcb3[count][4]=41;
        ++count;
        }

        if(p->flag3[IJKm1]<0)
        if(p->fgc[IJK][4]==0)
        {
        p->gcb3[count][0]=i;
        p->gcb3[count][1]=j;
        p->gcb3[count][2]=k;
        p->gcb3[count][3]=5;

        if(p->flag3[IJKm1]==TOPO)
        p->gcb3[count][4]=5;

        if(p->flag3[IJKm1]==OBJ)
        p->gcb3[count][4]=21;
        
        if(p->flag3[IJKm1]==SOLID)
        p->gcb3[count][4]=22;
		
		if(p->flag3[IJKm1]==FLT)
        p->gcb3[count][4]=41;
        ++count;
        }

        if(p->flag3[IJKp1]<0)
        if(p->fgc[IJK][5]==0)
        {
        p->gcb3[count][0]=i;
        p->gcb3[count][1]=j;
        p->gcb3[count][2]=k;
        p->gcb3[count][3]=6;

        if(p->flag3[IJKp1]==TOPO)
        p->gcb3[count][4]=5;

        if(p->flag3[IJKp1]==OBJ)
        p->gcb3[count][4]=21;
        
        if(p->flag3[IJKp1]==SOLID)
        p->gcb3[count][4]=22;

        if(p->flag3[IJKp1]==FLT)
        p->gcb3[count][4]=41;
        ++count;
        }
    }
    
    // extra parax + gcb
    for(q=0;q<p->gcpara1_count;++q)
    {
    i=p->gcpara1[q][0];
    j=p->gcpara1[q][1];
    k=p->gcpara1[q][2];
        
        if(p->flag3[IJK]==0 && p->flag3[Im1JK]==1)
        {
        p->gcb3[count][0]=i-1;
        p->gcb3[count][1]=j;
        p->gcb3[count][2]=k;
        p->gcb3[count][3]=4;

        if(p->flag3[IJK]==TOPO)
        p->gcb3[count][4]=5;

        if(p->flag3[IJK]==OBJ)
        p->gcb3[count][4]=21;
        
        if(p->flag3[IJK]==SOLID)
        p->gcb3[count][4]=22;
		
		if(p->flag3[IJK]==FLT)
        p->gcb3[count][4]=41;
        ++count;
        }
    }

    for(q=0;q<p->gcpara3_count;++q)
    {
    i=p->gcpara3[q][0];
    j=p->gcpara3[q][1];
    k=p->gcpara3[q][2];
        
        if(p->flag3[IJK]==0 && p->flag3[IJm1K]==1)
        {
        p->gcb3[count][0]=i;
        p->gcb3[count][1]=j-1;
        p->gcb3[count][2]=k;
        p->gcb3[count][3]=2;

        if(p->flag3[IJK]==TOPO)
        p->gcb3[count][4]=5;

        if(p->flag3[IJK]==OBJ)
        p->gcb3[count][4]=21;
        
        if(p->flag3[IJK]==SOLID)
        p->gcb3[count][4]=22;
		
		if(p->flag3[IJK]==FLT)
        p->gcb3[count][4]=41;
        ++count;
        }
    }

	for(q=0;q<p->gcpara5_count;++q)
	{
    i=p->gcpara5[q][0];
    j=p->gcpara5[q][1];
    k=p->gcpara5[q][2];
        
        if(p->flag3[IJK]==0 && p->flag3[IJKm1]==1)
        {
        p->gcb3[count][0]=i;
        p->gcb3[count][1]=j;
        p->gcb3[count][2]=k-1;
        p->gcb3[count][3]=6;

        if(p->flag3[IJK]==TOPO)
        p->gcb3[count][4]=5;

        if(p->flag3[IJK]==OBJ)
        p->gcb3[count][4]=21;
        
        if(p->flag3[IJK]==SOLID)
        p->gcb3[count][4]=22;
		
		if(p->flag3[IJK]==FLT)
        p->gcb3[count][4]=41;
        ++count;
        }
	}

	for(q=0;q<p->gcpara4_count;++q)
	{
    i=p->gcpara4[q][0]+1;
    j=p->gcpara4[q][1];
    k=p->gcpara4[q][2];
        
        if(p->flag3[IJK]==0 && p->flag3[Ip1JK]==1)
        {
        p->gcb3[count][0]=i+1;
        p->gcb3[count][1]=j;
        p->gcb3[count][2]=k;
        p->gcb3[count][3]=1;

        if(p->flag3[IJK]==TOPO)
        p->gcb3[count][4]=5;

        if(p->flag3[IJK]==OBJ)
        p->gcb3[count][4]=21;
        
        if(p->flag3[IJK]==SOLID)
        p->gcb3[count][4]=22;
		
		if(p->flag3[IJK]==FLT)
        p->gcb3[count][4]=41;
        ++count;
        }
	}

	for(q=0;q<p->gcpara2_count;++q)
	{
    i=p->gcpara2[q][0];
    j=p->gcpara2[q][1]+1;
    k=p->gcpara2[q][2];
        
        if(p->flag3[IJK]==0 && p->flag3[IJp1K]==1)
        {
        p->gcb3[count][0]=i;
        p->gcb3[count][1]=j+1;
        p->gcb3[count][2]=k;
        p->gcb3[count][3]=3;

        if(p->flag3[IJK]==TOPO)
        p->gcb3[count][4]=5;

        if(p->flag3[IJK]==OBJ)
        p->gcb3[count][4]=21;
        
        if(p->flag3[IJK]==SOLID)
        p->gcb3[count][4]=22;
		
		if(p->flag3[IJK]==FLT)
        p->gcb3[count][4]=41;
        ++count;
        }
	}

	for(q=0;q<p->gcpara6_count;++q)
	{
	i=p->gcpara6[q][0];
    j=p->gcpara6[q][1];
    k=p->gcpara6[q][2]+1;
        
        if(p->flag3[IJK]==0 && p->flag3[IJKp1]==1)
        {
        p->gcb3[count][0]=i;
        p->gcb3[count][1]=j;
        p->gcb3[count][2]=k+1;
        p->gcb3[count][3]=5;

        if(p->flag3[IJK]==TOPO)
        p->gcb3[count][4]=5;

        if(p->flag3[IJK]==OBJ)
        p->gcb3[count][4]=21;
        
        if(p->flag3[IJK]==SOLID)
        p->gcb3[count][4]=22;
		
		if(p->flag3[IJK]==FLT)
        p->gcb3[count][4]=41;
        ++count;
        }
	}
	
    for(q=p->gcb3_count;q<count;++q)
    {
		if(p->gcb3[q][3]==5 || p->gcb3[q][3]==6)
		p->gcd3[q]=0.5*p->DXM;

		if(p->gcb3[q][3]!=5 && p->gcb3[q][3]!=6)
		p->gcd3[q]=p->DXM;
    }

    p->gcb3_count=count;
	
	//cout<<p->mpirank<<" old gcb3: "<<p->gcb3_count<<" new gcb3: "<<count<<endl;


   // remove unactive gcb
	QGC3LOOP
	{
	    i=p->gcb3[q][0];
		j=p->gcb3[q][1];
		k=p->gcb3[q][2];
        if(p->flag3[IJK]<0)
        p->gcb3[q][3]=-fabs(p->gcb3[q][3]);
	}
    
    /* count1=0;
	QGCB3
	++count1;
    
    int count2=0;
	QGC3LOOP
	++count2;
    
    count=0;
    WLOOP
    {
        if(p->flag3[Im1JK]<0)
        ++count;
        
        if(p->flag3[Ip1JK]<0)
        ++count;
        
        if(p->flag3[IJm1K]<0)
        ++count;
        
        if(p->flag3[IJp1K]<0)
        ++count;
        
        if(p->flag3[IJKm1]<0)
        ++count;
        
        if(p->flag3[IJKp1]<0)
        ++count;
	}
	
	cout<<p->mpirank<<" GCB3: "<<p->gcb3_count<<" GCB3_direct_all: "<<count1<<" GCB3_direct: "<<count2<<" GCB3_LOOP: "<<count<<endl;
    */

	p->del_Iarray(p->fgc,imax*jmax*kmax,6);
}

