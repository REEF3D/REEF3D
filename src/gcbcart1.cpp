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

#include"mgc1.h"
#include"lexer.h"

void mgc1::fillgcb(lexer *p)
{
    int q,n;
	
	p->Iarray(p->fgc,imax*jmax*kmax,6);


	if(p->gcb1_count!=p->gcb4_count)
	{
	p->Iresize(p->gcb1,p->gcb1_count, p->gcb4_count, 6, 6); 	
	p->Dresize(p->gcd1,p->gcb1_count, p->gcb4_count); 
	
	p->gcb1_count=p->gcb4_count;
	}
	
	
	QGCB4
	{
	for(n=0;n<5;++n)
	p->gcb1[q][n]=p->gcb4[q][n];

    if(p->gcb1[q][3]==1 || p->gcb1[q][3]==4)
	p->gcd1[q]=p->gcd4[q];

	if(p->gcb1[q][3]!=1 && p->gcb1[q][3]!=4)
	p->gcd1[q]=p->gcd4[q];
	}
	

    QGC1LOOP
	{
	    i=p->gcb1[q][0];
		j=p->gcb1[q][1];
		k=p->gcb1[q][2];
        
        if(p->gcb1[q][3]==1 || p->gcb1[q][3]==4)
        {
        p->gcd1[q] += 0.5*p->DXP[IP];
        
        }

		p->fgc[IJK][p->gcb1[q][3]-1]=1;
	}


	QGC1LOOP
	{
	    i=p->gcb1[q][0];
		j=p->gcb1[q][1];
		k=p->gcb1[q][2];

            if(p->gcb1[q][3]==4 && (p->periodic1!=1 || i+p->origin_i<p->gknox-1))
            p->gcb1[q][0]-=1;
		
	}
    
    QGC1LOOP
	{
	    i=p->gcb1[q][0];
		j=p->gcb1[q][1];
		k=p->gcb1[q][2];

            if(p->gcb1[q][3]!=4 && p->fgc[IJK][3]==1 && (p->periodic1!=1 || i+p->origin_i<p->gknox-1))
            p->gcb1[q][3]=-fabs(p->gcb1[q][3]);
		
	}
}

void mgc1::extragcb(lexer *p)
{	
    int count;
    int q,n,l;

	for(n=0;n<imax*jmax*kmax;++n)
	for(q=0;q<6;++q)
	p->fgc[n][q]=0;

    QGC1LOOP
	{
	    i=p->gcb1[q][0];
		j=p->gcb1[q][1];
		k=p->gcb1[q][2];

		p->fgc[IJK][p->gcb1[q][3]-1]=1;
	}
	
    count=p->gcb1_count;
	
// Find number of extragcb

	ULOOP
    {
        if(p->flag1[UIm1JK]<0)
        if(p->fgc[IJK][0]==0)
        ++count;
		
		if(p->flag1[UIJp1K]<0)
        if(p->fgc[IJK][1]==0)
        ++count;

        if(p->flag1[UIJm1K]<0)
        if(p->fgc[IJK][2]==0)
        ++count;

        if(p->flag1[UIp1JK]<0)
        if(p->fgc[IJK][3]==0)
        ++count;

        if(p->flag1[UIJKm1]<0)
        if(p->fgc[IJK][4]==0)
        ++count;

        if(p->flag1[UIJKp1]<0)
        if(p->fgc[IJK][5]==0)
        ++count;
    }
    
    // extra parax + gcb
    for(q=0;q<p->gcpara1_count;++q)
    {
    i=p->gcpara1[q][0];
    j=p->gcpara1[q][1];
    k=p->gcpara1[q][2];
        
        if(p->flag1[UIJK]==0 && p->flag1[UIm1JK]==1)
        ++count;
    }

    for(q=0;q<p->gcpara3_count;++q)
    {
    i=p->gcpara3[q][0];
    j=p->gcpara3[q][1];
    k=p->gcpara3[q][2];
        
        if(p->flag1[UIJK]==0 && p->flag1[UIJm1K]==1)
        ++count;
    }

	for(q=0;q<p->gcpara5_count;++q)
	{
    i=p->gcpara5[q][0];
    j=p->gcpara5[q][1];
    k=p->gcpara5[q][2];
        
        if(p->flag1[UIJK]==0 && p->flag1[UIJKm1]==1)
        ++count;
	}

	for(q=0;q<p->gcpara4_count;++q)
	{
    i=p->gcpara4[q][0]+1;
    j=p->gcpara4[q][1];
    k=p->gcpara4[q][2];
        
        if(p->flag1[UIJK]==0 && p->flag1[UIp1JK]==1)
        ++count;
	}

	for(q=0;q<p->gcpara2_count;++q)
	{
    i=p->gcpara2[q][0];
    j=p->gcpara2[q][1]+1;
    k=p->gcpara2[q][2];
        
        if(p->flag1[UIJK]==0 && p->flag1[UIJp1K]==1)
        ++count;
	}

	for(q=0;q<p->gcpara6_count;++q)
	{
	i=p->gcpara6[q][0];
    j=p->gcpara6[q][1];
    k=p->gcpara6[q][2]+1;
        
        if(p->flag1[UIJK]==0 && p->flag1[UIJKp1]==1)
        ++count;
	}
	//cout<<p->mpirank<<" old gcb1: "<<p->gcb1_count<<" new gcb1: "<<count<<endl;
	
	if(p->gcb1_count!=count)
	{
	p->Iresize(p->gcb1,p->gcb1_count, count, 6, 6); 
	p->Dresize(p->gcd1,p->gcb1_count, count); 
	}

// Store extragcb in gcb vec
	
	count=p->gcb1_count;
	
    ULOOP
    {	
        if(p->flag1[UIm1JK]<0)
        if(p->fgc[IJK][0]==0)
        {
        p->gcb1[count][0]=i;
        p->gcb1[count][1]=j;
        p->gcb1[count][2]=k;
        p->gcb1[count][3]=1;

        if(p->flag1[UIm1JK]==TOPO)
        p->gcb1[count][4]=5;

        if(p->flag1[UIm1JK]==OBJ)
        p->gcb1[count][4]=21;
        
        if(p->flag1[UIm1JK]==SOLID)
        p->gcb1[count][4]=22;
		
		if(p->flag1[UIm1JK]==FLT)
        p->gcb1[count][4]=41;
        ++count;
        }
		
		if(p->flag1[UIJp1K]<0)
        if(p->fgc[IJK][1]==0)
        {
        p->gcb1[count][0]=i;
        p->gcb1[count][1]=j;
        p->gcb1[count][2]=k;
        p->gcb1[count][3]=2;

        if(p->flag1[UIJp1K]==TOPO)
        p->gcb1[count][4]=5;

        if(p->flag1[UIJp1K]==OBJ)
        p->gcb1[count][4]=21;
        
        if(p->flag1[UIJp1K]==SOLID)
        p->gcb1[count][4]=22;
		
		if(p->flag1[UIJp1K]==FLT)
        p->gcb1[count][4]=41;
        ++count;
        }

        if(p->flag1[UIJm1K]<0)
        if(p->fgc[IJK][2]==0)
        {
        p->gcb1[count][0]=i;
        p->gcb1[count][1]=j;
        p->gcb1[count][2]=k;
        p->gcb1[count][3]=3;

        if(p->flag1[UIJm1K]==TOPO)
        p->gcb1[count][4]=5;

        if(p->flag1[UIJm1K]==OBJ)
        p->gcb1[count][4]=21;
        
        if(p->flag1[UIJm1K]==SOLID)
        p->gcb1[count][4]=22;
		
		if(p->flag1[UIJm1K]==FLT)
        p->gcb1[count][4]=41;
        ++count;
        }

        if(p->flag1[UIp1JK]<0)
        if(p->fgc[IJK][3]==0)
        {
        p->gcb1[count][0]=i;
        p->gcb1[count][1]=j;
        p->gcb1[count][2]=k;
        p->gcb1[count][3]=4;

        if(p->flag1[UIp1JK]==TOPO)
        p->gcb1[count][4]=5;

        if(p->flag1[UIp1JK]==OBJ)
        p->gcb1[count][4]=21;
        
        if(p->flag1[UIp1JK]==SOLID)
        p->gcb1[count][4]=22;
		
		if(p->flag1[UIp1JK]==FLT)
        p->gcb1[count][4]=41;
        ++count;
        }

        if(p->flag1[UIJKm1]<0)
        if(p->fgc[IJK][4]==0)
        {
        p->gcb1[count][0]=i;
        p->gcb1[count][1]=j;
        p->gcb1[count][2]=k;
        p->gcb1[count][3]=5;

        if(p->flag1[UIJKm1]==TOPO)
        p->gcb1[count][4]=5;

        if(p->flag1[UIJKm1]==OBJ)
        p->gcb1[count][4]=21;
        
        if(p->flag1[UIJKm1]==SOLID)
        p->gcb1[count][4]=22;
		
		if(p->flag1[UIJKm1]==FLT)
        p->gcb1[count][4]=41;

        ++count;
        }

        if(p->flag1[UIJKp1]<0)
        if(p->fgc[IJK][5]==0)
        {
        p->gcb1[count][0]=i;
        p->gcb1[count][1]=j;
        p->gcb1[count][2]=k;
        p->gcb1[count][3]=6;

        if(p->flag1[UIJKp1]==TOPO)
        p->gcb1[count][4]=5;

        if(p->flag1[UIJKp1]==OBJ)
        p->gcb1[count][4]=21;
        
        if(p->flag1[UIJKp1]==SOLID)
        p->gcb1[count][4]=22;
		
		if(p->flag1[UIJKp1]==FLT)
        p->gcb1[count][4]=41;
        ++count;
        }
    }
    
    
    // extra parax + gcb
    for(q=0;q<p->gcpara1_count;++q)
    {
    i=p->gcpara1[q][0];
    j=p->gcpara1[q][1];
    k=p->gcpara1[q][2];
        
        if(p->flag1[UIJK]==0 && p->flag1[UIm1JK]==1)
        {
        p->gcb1[count][0]=i-1;
        p->gcb1[count][1]=j;
        p->gcb1[count][2]=k;
        p->gcb1[count][3]=4;

        if(p->flag1[UIJK]==TOPO)
        p->gcb1[count][4]=5;

        if(p->flag1[UIJK]==OBJ)
        p->gcb1[count][4]=21;
        
        if(p->flag1[UIJK]==SOLID)
        p->gcb1[count][4]=22;
		
		if(p->flag1[UIJK]==FLT)
        p->gcb1[count][4]=41;
        ++count;
        }
    }

    for(q=0;q<p->gcpara3_count;++q)
    {
    i=p->gcpara3[q][0];
    j=p->gcpara3[q][1];
    k=p->gcpara3[q][2];
        
        if(p->flag1[UIJK]==0 && p->flag1[UIJm1K]==1)
        {
        p->gcb1[count][0]=i;
        p->gcb1[count][1]=j-1;
        p->gcb1[count][2]=k;
        p->gcb1[count][3]=2;

        if(p->flag1[UIJK]==TOPO)
        p->gcb1[count][4]=5;

        if(p->flag1[UIJK]==OBJ)
        p->gcb1[count][4]=21;
        
        if(p->flag1[UIJK]==SOLID)
        p->gcb1[count][4]=22;
		
		if(p->flag1[UIJK]==FLT)
        p->gcb1[count][4]=41;
        ++count;
        }
    }

	for(q=0;q<p->gcpara5_count;++q)
	{
    i=p->gcpara5[q][0];
    j=p->gcpara5[q][1];
    k=p->gcpara5[q][2];
        
        if(p->flag1[UIJK]==0 && p->flag1[UIJKm1]==1)
        {
        p->gcb1[count][0]=i;
        p->gcb1[count][1]=j;
        p->gcb1[count][2]=k-1;
        p->gcb1[count][3]=6;

        if(p->flag1[UIJK]==TOPO)
        p->gcb1[count][4]=5;

        if(p->flag1[UIJK]==OBJ)
        p->gcb1[count][4]=21;
        
        if(p->flag1[UIJK]==SOLID)
        p->gcb1[count][4]=22;
		
		if(p->flag1[UIJK]==FLT)
        p->gcb1[count][4]=41;
        ++count;
        }
	}

	for(q=0;q<p->gcpara4_count;++q)
	{
    i=p->gcpara4[q][0]+1;
    j=p->gcpara4[q][1];
    k=p->gcpara4[q][2];
        
        if(p->flag1[UIJK]==0 && p->flag1[UIp1JK]==1)
        {
        p->gcb1[count][0]=i+1;
        p->gcb1[count][1]=j;
        p->gcb1[count][2]=k;
        p->gcb1[count][3]=1;

        if(p->flag1[UIJK]==TOPO)
        p->gcb1[count][4]=5;

        if(p->flag1[UIJK]==OBJ)
        p->gcb1[count][4]=21;
        
        if(p->flag1[UIJK]==SOLID)
        p->gcb1[count][4]=22;
		
		if(p->flag1[UIJK]==FLT)
        p->gcb1[count][4]=41;
        ++count;
        }
	}

	for(q=0;q<p->gcpara2_count;++q)
	{
    i=p->gcpara2[q][0];
    j=p->gcpara2[q][1]+1;
    k=p->gcpara2[q][2];
        
        if(p->flag1[UIJK]==0 && p->flag1[UIJp1K]==1)
        {
        p->gcb1[count][0]=i;
        p->gcb1[count][1]=j+1;
        p->gcb1[count][2]=k;
        p->gcb1[count][3]=3;

        if(p->flag1[UIJK]==TOPO)
        p->gcb1[count][4]=5;

        if(p->flag1[UIJK]==OBJ)
        p->gcb1[count][4]=21;
        
        if(p->flag1[UIJK]==SOLID)
        p->gcb1[count][4]=22;
		
		if(p->flag1[UIJK]==FLT)
        p->gcb1[count][4]=41;
        ++count;
        }
	}

	for(q=0;q<p->gcpara6_count;++q)
	{
	i=p->gcpara6[q][0];
    j=p->gcpara6[q][1];
    k=p->gcpara6[q][2]+1;
        
        if(p->flag1[UIJK]==0 && p->flag1[UIJKp1]==1)
        {
        p->gcb1[count][0]=i;
        p->gcb1[count][1]=j;
        p->gcb1[count][2]=k+1;
        p->gcb1[count][3]=5;

        if(p->flag1[UIJK]==TOPO)
        p->gcb1[count][4]=5;

        if(p->flag1[UIJK]==OBJ)
        p->gcb1[count][4]=21;
        
        if(p->flag1[UIJK]==SOLID)
        p->gcb1[count][4]=22;
		
		if(p->flag1[UIJK]==FLT)
        p->gcb1[count][4]=41;
        ++count;
        }
	}
	

    for(q=p->gcb1_count;q<count;++q)
    {
		if(p->gcb1[q][3]==1 || p->gcb1[q][3]==4)
		p->gcd1[q]=0.5*p->DXM;

		if(p->gcb1[q][3]!=1 && p->gcb1[q][3]!=4)
		p->gcd1[q]=p->DXM;
    }
	
    p->gcb1_count=count;
	
   
    // remove unactive gcb
	QGC1LOOP
	{
	    i=p->gcb1[q][0];
		j=p->gcb1[q][1];
		k=p->gcb1[q][2];
        if(p->flag1[UIJK]<0)
        p->gcb1[q][3]=-fabs(p->gcb1[q][3]);
	}
    
	p->del_Iarray(p->fgc,imax*jmax*kmax,6);
	
}


