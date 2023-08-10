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

#include"mgc2.h"
#include"lexer.h"

void mgc2::fillgcb(lexer *p)
{
    int q,n;
	
	p->Iarray(p->fgc,imax*jmax*kmax,6);
	
// ----

	if(p->gcb2_count!=p->gcb4_count)
	{
	p->Iresize(p->gcb2,p->gcb2_count, p->gcb4_count, 6, 6); 
	p->Dresize(p->gcd2,p->gcb2_count, p->gcb4_count); 
	
	p->gcb2_count=p->gcb4_count;
	}

	QGCB4
	{
	for(n=0;n<5;++n)
	p->gcb2[q][n]=p->gcb4[q][n];

	if(p->gcb2[q][3]==2 || p->gcb2[q][3]==3)
	p->gcd2[q]=p->gcd4[q];

	if(p->gcb2[q][3]!=2 && p->gcb2[q][3]!=3)
	p->gcd2[q]=p->gcd4[q];
	}

    QGC2LOOP
	{
	    i=p->gcb2[q][0];
		j=p->gcb2[q][1];
		k=p->gcb2[q][2];
        
        if(p->gcb2[q][3]==2 || p->gcb2[q][3]==3)
        p->gcd2[q] += 0.5*p->DYP[JP];

		p->fgc[IJK][p->gcb2[q][3]-1]=1;
	}
	
	
	
	QGC2LOOP
	{
	    i=p->gcb2[q][0];
		j=p->gcb2[q][1];
		k=p->gcb2[q][2];

            if(p->gcb2[q][3]==2 && (p->periodic2!=1 || j+p->origin_j<p->gknoy-1))
            p->gcb2[q][1]-=1;
	}
    
    QGC2LOOP
	{
	    i=p->gcb2[q][0];
		j=p->gcb2[q][1];
		k=p->gcb2[q][2];

            if(p->gcb2[q][3]!=2 && p->fgc[IJK][1]==1 && (p->periodic2!=1 || j+p->origin_j<p->gknoy-1))
            p->gcb2[q][3]=-fabs(p->gcb2[q][3]);
	}
}

void mgc2::extragcb(lexer *p)
{
    int count;
    int q,n;

	for(n=0;n<imax*jmax*kmax;++n)
	for(int l=0;l<6;++l)
	p->fgc[n][l]=0;

    QGC2LOOP
	{
	    i=p->gcb2[q][0];
		j=p->gcb2[q][1];
		k=p->gcb2[q][2];

		p->fgc[IJK][p->gcb2[q][3]-1]=1;
	}

    count=p->gcb2_count;
	
// Find number of extragcb
	
	VLOOP
    {
        if(p->flag2[Im1JK]<0)
        if(p->fgc[IJK][0]==0)
        ++count;

        if(p->flag2[IJm1K]<0)
        if(p->fgc[IJK][2]==0)
        ++count;

        if(p->flag2[IJp1K]<0)
        if(p->fgc[IJK][1]==0)
        ++count;

        if(p->flag2[Ip1JK]<0)
        if(p->fgc[IJK][3]==0)
        ++count;

        if(p->flag2[IJKm1]<0)
        if(p->fgc[IJK][4]==0)
        ++count;

        if(p->flag2[IJKp1]<0)
        if(p->fgc[IJK][5]==0)
        ++count;
    }
    
    // extra parax + gcb
    for(q=0;q<p->gcpara1_count;++q)
    {
    i=p->gcpara1[q][0];
    j=p->gcpara1[q][1];
    k=p->gcpara1[q][2];
        
        if(p->flag2[IJK]==0 && p->flag2[Im1JK]==1)
        ++count;
    }

    for(q=0;q<p->gcpara3_count;++q)
    {
    i=p->gcpara3[q][0];
    j=p->gcpara3[q][1];
    k=p->gcpara3[q][2];
        
        if(p->flag2[IJK]==0 && p->flag2[IJm1K]==1)
        ++count;
    }

	for(q=0;q<p->gcpara5_count;++q)
	{
    i=p->gcpara5[q][0];
    j=p->gcpara5[q][1];
    k=p->gcpara5[q][2];
        
        if(p->flag2[IJK]==0 && p->flag2[IJKm1]==1)
        ++count;
	}

	for(q=0;q<p->gcpara4_count;++q)
	{
    i=p->gcpara4[q][0]+1;
    j=p->gcpara4[q][1];
    k=p->gcpara4[q][2];
        
        if(p->flag2[IJK]==0 && p->flag2[Ip1JK]==1)
        ++count;
	}

	for(q=0;q<p->gcpara2_count;++q)
	{
    i=p->gcpara2[q][0];
    j=p->gcpara2[q][1]+1;
    k=p->gcpara2[q][2];
        
        if(p->flag2[IJK]==0 && p->flag2[IJp1K]==1)
        ++count;
	}

	for(q=0;q<p->gcpara6_count;++q)
	{
	i=p->gcpara6[q][0];
    j=p->gcpara6[q][1];
    k=p->gcpara6[q][2]+1;
        
        if(p->flag2[IJK]==0 && p->flag2[IJKp1]==1)
        ++count;
	}

	if(p->gcb2_count!=count)
	{
	p->Iresize(p->gcb2,p->gcb2_count, count, 6, 6); 
	p->Dresize(p->gcd2,p->gcb2_count, count); 
	}


// Store extragcb in gcb vec
	
	count=p->gcb2_count;
	
    VLOOP
    {
        if(p->flag2[Im1JK]<0)
        if(p->fgc[IJK][0]==0)
        {
        p->gcb2[count][0]=i;
        p->gcb2[count][1]=j;
        p->gcb2[count][2]=k;
        p->gcb2[count][3]=1;

        if(p->flag2[Im1JK]==TOPO)
        p->gcb2[count][4]=5;

        if(p->flag2[Im1JK]==OBJ)
        p->gcb2[count][4]=21;
        
        if(p->flag2[Im1JK]==SOLID)
        p->gcb2[count][4]=22;
        
        if(p->flag2[Im1JK]==FLT)
        p->gcb2[count][4]=41;
        ++count;
        }

        if(p->flag2[IJm1K]<0)
        if(p->fgc[IJK][2]==0)
        {
        p->gcb2[count][0]=i;
        p->gcb2[count][1]=j;
        p->gcb2[count][2]=k;
        p->gcb2[count][3]=3;

        if(p->flag2[IJm1K]==TOPO)
        p->gcb2[count][4]=5;

        if(p->flag2[IJm1K]==OBJ)
        p->gcb2[count][4]=21;
        
        if(p->flag2[IJm1K]==SOLID)
        p->gcb2[count][4]=22;
        
        if(p->flag2[IJm1K]==FLT)
        p->gcb2[count][4]=41;    
        ++count;
        }

        if(p->flag2[IJp1K]<0)
        if(p->fgc[IJK][1]==0)
        {
        p->gcb2[count][0]=i;
        p->gcb2[count][1]=j;
        p->gcb2[count][2]=k;
        p->gcb2[count][3]=2;

        if(p->flag2[IJp1K]==TOPO)
        p->gcb2[count][4]=5;

        if(p->flag2[IJp1K]==OBJ)
        p->gcb2[count][4]=21;
        
        if(p->flag2[IJp1K]==SOLID)
        p->gcb2[count][4]=22;
        
        if(p->flag2[IJp1K]==FLT)
        p->gcb2[count][4]=41;
        ++count;
        }

        if(p->flag2[Ip1JK]<0)
        if(p->fgc[IJK][3]==0)
        {
        p->gcb2[count][0]=i;
        p->gcb2[count][1]=j;
        p->gcb2[count][2]=k;
        p->gcb2[count][3]=4;

        if(p->flag2[Ip1JK]==TOPO)
        p->gcb2[count][4]=5;

        if(p->flag2[Ip1JK]==OBJ)
        p->gcb2[count][4]=21;
        
        if(p->flag2[Ip1JK]==SOLID)
        p->gcb2[count][4]=22;
        
        if(p->flag2[Ip1JK]==FLT)
        p->gcb2[count][4]=41;
        ++count;
        }

        if(p->flag2[IJKm1]<0)
        if(p->fgc[IJK][4]==0)
        {
        p->gcb2[count][0]=i;
        p->gcb2[count][1]=j;
        p->gcb2[count][2]=k;
        p->gcb2[count][3]=5;

        if(p->flag2[IJKm1]==TOPO)
        p->gcb2[count][4]=5;

        if(p->flag2[IJKm1]==OBJ)
        p->gcb2[count][4]=21;
        
        if(p->flag2[IJKm1]==SOLID)
        p->gcb2[count][4]=22;

        if(p->flag2[IJKm1]==FLT)
        p->gcb2[count][4]=41;
        ++count;
        }

        if(p->flag2[IJKp1]<0)
        if(p->fgc[IJK][5]==0)
        {
        p->gcb2[count][0]=i;
        p->gcb2[count][1]=j;
        p->gcb2[count][2]=k;
        p->gcb2[count][3]=6;

        if(p->flag2[IJKp1]==TOPO)
        p->gcb2[count][4]=5;

        if(p->flag2[IJKp1]==OBJ)
        p->gcb2[count][4]=21;
        
        if(p->flag2[IJKp1]==SOLID)
        p->gcb2[count][4]=22;
        
        if(p->flag2[IJKp1]==FLT)
        p->gcb2[count][4]=41;
        ++count;
        }
    }
    
    // extra parax + gcb
    for(q=0;q<p->gcpara1_count;++q)
    {
    i=p->gcpara1[q][0];
    j=p->gcpara1[q][1];
    k=p->gcpara1[q][2];
        
        if(p->flag2[IJK]==0 && p->flag2[Im1JK]==1)
        {
        p->gcb2[count][0]=i-1;
        p->gcb2[count][1]=j;
        p->gcb2[count][2]=k;
        p->gcb2[count][3]=4;

        if(p->flag2[IJK]==TOPO)
        p->gcb2[count][4]=5;

        if(p->flag2[IJK]==OBJ)
        p->gcb2[count][4]=21;
        
        if(p->flag2[IJK]==SOLID)
        p->gcb2[count][4]=22;
		
		if(p->flag2[IJK]==FLT)
        p->gcb2[count][4]=41;
        ++count;
        }
    }

    for(q=0;q<p->gcpara3_count;++q)
    {
    i=p->gcpara3[q][0];
    j=p->gcpara3[q][1];
    k=p->gcpara3[q][2];
        
        if(p->flag2[IJK]==0 && p->flag2[IJm1K]==1)
        {
        p->gcb2[count][0]=i;
        p->gcb2[count][1]=j-1;
        p->gcb2[count][2]=k;
        p->gcb2[count][3]=2;

        if(p->flag2[IJK]==TOPO)
        p->gcb2[count][4]=5;

        if(p->flag2[IJK]==OBJ)
        p->gcb2[count][4]=21;
        
        if(p->flag2[IJK]==SOLID)
        p->gcb2[count][4]=22;
		
		if(p->flag2[IJK]==FLT)
        p->gcb2[count][4]=41;
        ++count;
        }
    }

	for(q=0;q<p->gcpara5_count;++q)
	{
    i=p->gcpara5[q][0];
    j=p->gcpara5[q][1];
    k=p->gcpara5[q][2];
        
        if(p->flag2[IJK]==0 && p->flag2[IJKm1]==1)
        {
        p->gcb2[count][0]=i;
        p->gcb2[count][1]=j;
        p->gcb2[count][2]=k-1;
        p->gcb2[count][3]=6;

        if(p->flag2[IJK]==TOPO)
        p->gcb2[count][4]=5;

        if(p->flag2[IJK]==OBJ)
        p->gcb2[count][4]=21;
        
        if(p->flag2[IJK]==SOLID)
        p->gcb2[count][4]=22;
		
		if(p->flag2[IJK]==FLT)
        p->gcb2[count][4]=41;
        ++count;
        }
	}

	for(q=0;q<p->gcpara4_count;++q)
	{
    i=p->gcpara4[q][0]+1;
    j=p->gcpara4[q][1];
    k=p->gcpara4[q][2];
        
        if(p->flag2[IJK]==0 && p->flag2[Ip1JK]==1)
        {
        p->gcb2[count][0]=i+1;
        p->gcb2[count][1]=j;
        p->gcb2[count][2]=k;
        p->gcb2[count][3]=1;

        if(p->flag2[IJK]==TOPO)
        p->gcb2[count][4]=5;

        if(p->flag2[IJK]==OBJ)
        p->gcb2[count][4]=21;
        
        if(p->flag2[IJK]==SOLID)
        p->gcb2[count][4]=22;
		
		if(p->flag2[IJK]==FLT)
        p->gcb2[count][4]=41;
        ++count;
        }
	}

	for(q=0;q<p->gcpara2_count;++q)
	{
    i=p->gcpara2[q][0];
    j=p->gcpara2[q][1]+1;
    k=p->gcpara2[q][2];
        
        if(p->flag2[IJK]==0 && p->flag2[IJp1K]==1)
        {
        p->gcb2[count][0]=i;
        p->gcb2[count][1]=j+1;
        p->gcb2[count][2]=k;
        p->gcb2[count][3]=3;

        if(p->flag2[IJK]==TOPO)
        p->gcb2[count][4]=5;

        if(p->flag2[IJK]==OBJ)
        p->gcb2[count][4]=21;
        
        if(p->flag2[IJK]==SOLID)
        p->gcb2[count][4]=22;
		
		if(p->flag2[IJK]==FLT)
        p->gcb2[count][4]=41;
        ++count;
        }
	}

	for(q=0;q<p->gcpara6_count;++q)
	{
	i=p->gcpara6[q][0];
    j=p->gcpara6[q][1];
    k=p->gcpara6[q][2]+1;
        
        if(p->flag2[IJK]==0 && p->flag2[IJKp1]==1)
        {
        p->gcb2[count][0]=i;
        p->gcb2[count][1]=j;
        p->gcb2[count][2]=k+1;
        p->gcb2[count][3]=5;

        if(p->flag2[IJK]==TOPO)
        p->gcb2[count][4]=5;

        if(p->flag2[IJK]==OBJ)
        p->gcb2[count][4]=21;
        
        if(p->flag2[IJK]==SOLID)
        p->gcb2[count][4]=22;
		
		if(p->flag2[IJK]==FLT)
        p->gcb2[count][4]=41;
        ++count;
        }
	}
	
    for(q=p->gcb2_count;q<count;++q)
    {
		if(p->gcb2[q][3]==2 || p->gcb2[q][3]==3)
		p->gcd2[q]=0.5*p->DXM;

		if(p->gcb2[q][3]!=2 && p->gcb2[q][3]!=3)
		p->gcd2[q]=p->DXM;
    }
	
    p->gcb2_count=count;

    // remove unactive gcb
	QGC2LOOP
	{
	    i=p->gcb2[q][0];
		j=p->gcb2[q][1];
		k=p->gcb2[q][2];
        if(p->flag2[IJK]<0)
        p->gcb2[q][3]=-fabs(p->gcb2[q][3]);
	}
	
    /*
    int count1=0;
	QGCB2
	++count1;
    
    int count2=0;
	QGC2LOOP
	++count2;
    
    count=0;
    VLOOP
    {
        if(p->flag2[Im1JK]<0)
        ++count;
        
        if(p->flag2[Ip1JK]<0)
        ++count;
        
        if(p->flag2[IJm1K]<0)
        ++count;
        
        if(p->flag2[IJp1K]<0)
        ++count;
        
        if(p->flag2[IJKm1]<0)
        ++count;
        
        if(p->flag2[IJKp1]<0)
        ++count;
	}
	
	cout<<p->mpirank<<" GCB2: "<<p->gcb2_count<<" GCB2_direct_all: "<<count1<<" GCB2_direct: "<<count2<<" GCB2_LOOP: "<<count<<endl;
    */
	
	p->del_Iarray(p->fgc,imax*jmax*kmax,6);
}

