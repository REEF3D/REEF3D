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

#include"mgc4.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fieldint5.h"

mgc4::mgc4(lexer *p)
{
	imin=p->imin;
    imax=p->imax;
    jmin=p->jmin;
    jmax=p->jmax;
    kmin=p->kmin;
    kmax=p->kmax;
}

mgc4::~mgc4()
{
}

void mgc4::makemgc(lexer* p)
{
	p->Iarray(p->mgc4,kmax*jmax*imax);

//make gcdir
	p->gcdirsize4=1;	
	p->Iarray(p->gcorig4, p->gcdirsize4, 6,4);
}

void mgc4::mgcsetup(lexer* p)
{
	for(i=0;i<imax*jmax*kmax;++i)
	p->mgc4[i]=0;

	LOOP
	p->mgc4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin]=1;
}

void mgc4::fillmgc(lexer* p)
{
	int q,n;
	
//--------------------------
//WALL1
	QGC4LOOP
	{
        i=p->gcb4[q][0];
        j=p->gcb4[q][1];
        k=p->gcb4[q][2];
		
		if(p->gcb4[q][3]==1)
		for(n=0;n<p->margin;++n)
        p->mgc4[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]+=1;

		if(p->gcb4[q][3]==4)
		for(n=0;n<p->margin;++n)
		p->mgc4[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]+=1;

		if(p->gcb4[q][3]==3)
		for(n=0;n<p->margin;++n)
        p->mgc4[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]+=1;

		if(p->gcb4[q][3]==2)
		for(n=0;n<p->margin;++n)
		p->mgc4[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]+=1;

		if(p->gcb4[q][3]==5)
		for(n=0;n<p->margin;++n)
		p->mgc4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]+=1;

		if(p->gcb4[q][3]==6)
		for(n=0;n<p->margin;++n)
		p->mgc4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]+=1;
	}
    
    
//PARA1
    for(q=0;q<p->gcpara1_count;++q)
    {
    i=p->gcpara1[q][0];
    j=p->gcpara1[q][1];
    k=p->gcpara1[q][2];
        
        if(p->gcpara1[q][6]==1)
        for(n=0;n<p->margin;++n)
        p->mgc4[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]+=1;
    }
    
    for(q=0;q<p->gcpara4_count;++q)
	{
    i=p->gcpara4[q][0];
    j=p->gcpara4[q][1];
    k=p->gcpara4[q][2];
        
        if(p->gcpara4[q][6]==1)
        for(n=0;n<p->margin;++n)
		p->mgc4[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]+=1;
	}

    for(q=0;q<p->gcpara3_count;++q)
    {
    i=p->gcpara3[q][0];
    j=p->gcpara3[q][1];
    k=p->gcpara3[q][2];
        
        if(p->gcpara3[q][6]==1)
        for(n=0;n<p->margin;++n)
        p->mgc4[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]+=1;
    }
    
    for(q=0;q<p->gcpara2_count;++q)
	{
    i=p->gcpara2[q][0];
    j=p->gcpara2[q][1];
    k=p->gcpara2[q][2];
        
        if(p->gcpara2[q][6]==1)
        for(n=0;n<p->margin;++n)
        p->mgc4[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]+=1;
	}

	for(q=0;q<p->gcpara5_count;++q)
	{
    i=p->gcpara5[q][0];
    j=p->gcpara5[q][1];
    k=p->gcpara5[q][2];
        
        if(p->gcpara5[q][6]==1)
        for(n=0;n<p->margin;++n)
        p->mgc4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]+=1;
	}

	for(q=0;q<p->gcpara6_count;++q)
	{
	i=p->gcpara6[q][0];
    j=p->gcpara6[q][1];
    k=p->gcpara6[q][2];
        
        if(p->gcpara6[q][6]==1)
        for(n=0;n<p->margin;++n)
        p->mgc4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]+=1;
	}

//--------------------------
//WALL2
    p->gcextra4=10;
    
	QGC4LOOP
	{
        i=p->gcb4[q][0];
        j=p->gcb4[q][1];
        k=p->gcb4[q][2];

		if(p->gcb4[q][3]==1)
		for(n=0;n<p->margin;++n)
		if(p->mgc4[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]>1
		&& p->mgc4[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]<10)
        {
			p->mgc4[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]=p->gcextra4;
			++p->gcextra4;
        }

		if(p->gcb4[q][3]==4)
		for(n=0;n<p->margin;++n)
        if(p->mgc4[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]>1
		&& p->mgc4[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]<10)
        {
			p->mgc4[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]=p->gcextra4;
			++p->gcextra4;
        }

		if(p->gcb4[q][3]==3)
		for(n=0;n<p->margin;++n)
		if(p->mgc4[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]>1
		&& p->mgc4[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]<10)
        {
			p->mgc4[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]=p->gcextra4;
			++p->gcextra4;
        }

		if(p->gcb4[q][3]==2)
		for(n=0;n<p->margin;++n)
		if(p->mgc4[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]>1
		&& p->mgc4[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]<10)
        {
			p->mgc4[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]=p->gcextra4;
			++p->gcextra4;
        }

		if(p->gcb4[q][3]==5)
		for(n=0;n<p->margin;++n)
		if( p->mgc4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]>1
		&& p->mgc4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]<10)
        {
			p->mgc4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]=p->gcextra4;
			++p->gcextra4;
        }

		if(p->gcb4[q][3]==6)
		for(n=0;n<p->margin;++n)
		if(p->mgc4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]>1
		&& p->mgc4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]<10)
        {
			p->mgc4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]=p->gcextra4;
			++p->gcextra4;
        }
	}

//PARA2
    for(q=0;q<p->gcpara1_count;++q)
    {
    i=p->gcpara1[q][0];
    j=p->gcpara1[q][1];
    k=p->gcpara1[q][2];
        
        if(p->gcpara1[q][6]==1)
        for(n=0;n<p->margin;++n)
        if(p->mgc4[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]>1
		&& p->mgc4[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]<10)
        {
			p->mgc4[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]=p->gcextra4;
			++p->gcextra4;
        }
    }
    
    for(q=0;q<p->gcpara4_count;++q)
	{
    i=p->gcpara4[q][0];
    j=p->gcpara4[q][1];
    k=p->gcpara4[q][2];
        
        if(p->gcpara4[q][6]==1)
        for(n=0;n<p->margin;++n)
		if(p->mgc4[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]>1
		&& p->mgc4[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]<10)
        {
			p->mgc4[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]=p->gcextra4;
			++p->gcextra4;
        }
	}

    for(q=0;q<p->gcpara3_count;++q)
    {
    i=p->gcpara3[q][0];
    j=p->gcpara3[q][1];
    k=p->gcpara3[q][2];
        
        if(p->gcpara3[q][6]==1)
        for(n=0;n<p->margin;++n)
        if(p->mgc4[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]>1
		&& p->mgc4[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]<10)
        {
			p->mgc4[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]=p->gcextra4;
			++p->gcextra4;
        }
    }
    
    for(q=0;q<p->gcpara2_count;++q)
	{
    i=p->gcpara2[q][0];
    j=p->gcpara2[q][1];
    k=p->gcpara2[q][2];
        
        if(p->gcpara2[q][6]==1)
        for(n=0;n<p->margin;++n)
        if(p->mgc4[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]>1
		&& p->mgc4[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]<10)
        {
			p->mgc4[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]=p->gcextra4;
			++p->gcextra4;
        }
	}

	for(q=0;q<p->gcpara5_count;++q)
	{
    i=p->gcpara5[q][0];
    j=p->gcpara5[q][1];
    k=p->gcpara5[q][2];
        
        if(p->gcpara5[q][6]==1)
        for(n=0;n<p->margin;++n)
        if( p->mgc4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]>1
		&& p->mgc4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]<10)
        {
			p->mgc4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]=p->gcextra4;
			++p->gcextra4;
        }
	}

	for(q=0;q<p->gcpara6_count;++q)
	{
	i=p->gcpara6[q][0];
    j=p->gcpara6[q][1];
    k=p->gcpara6[q][2];
        
        if(p->gcpara6[q][6]==1)
        for(n=0;n<p->margin;++n)
        if(p->mgc4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]>1
		&& p->mgc4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]<10)
        {
			p->mgc4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]=p->gcextra4;
			++p->gcextra4;
        }
	}
    
    //cout<<p->mpirank<<" GCEXTRA4: "<<p->gcextra4<<" solid_est_gcextra: "<<MAX(MAX(p->solid_gcbextra_est,p->topo_gcbextra_est),p->tot_gcbextra_est)+10 <<endl;
}

void mgc4::gcdirfill(lexer* p)
{
// GCORIG
    int q,n;
    
	p->Iresize(p->gcorig4,p->gcdirsize4, p->gcextra4, 6, 6, 4, 4); 
	p->gcdirsize4=p->gcextra4;
	
	
	for(n=0;n<p->gcdirsize4;++n)
	for(q=0;q<6;++q)	
	for(qn=0;qn<4;++qn)	
	p->gcorig4[n][q][qn]=0;	

// WALL fill	
	QGC4LOOP
	{
        i=p->gcb4[q][0];
        j=p->gcb4[q][1];
        k=p->gcb4[q][2];

		if(p->gcb4[q][3]==1)
		for(n=0;n<p->margin;++n)
		if( p->mgc4[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]>1)
        {
			di = (n+1);
			p->gcorig4[p->mgc4[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]-10][0][di]=1;
        }

		if(p->gcb4[q][3]==4)
		for(n=0;n<p->margin;++n)
        if( p->mgc4[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]>1)
        {
			di = (n+1);
			p->gcorig4[p->mgc4[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]-10][3][di]=1;
        }

		if(p->gcb4[q][3]==3)
		for(n=0;n<p->margin;++n)
		if(p->mgc4[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]>1)
        {
			dj = (n+1);
			p->gcorig4[p->mgc4[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]-10][2][dj]=1;
	    }

		if(p->gcb4[q][3]==2)
		for(n=0;n<p->margin;++n)
		if( p->mgc4[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]>1)
        {
			dj = (n+1);
			p->gcorig4[p->mgc4[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]-10][1][dj]=1;
        }

		if(p->gcb4[q][3]==5)
		for(n=0;n<p->margin;++n)
		if( p->mgc4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]>1)
        {
			dk = (n+1);
			p->gcorig4[p->mgc4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]-10][4][dk]=1;
        }

		if(p->gcb4[q][3]==6)
		for(n=0;n<p->margin;++n)
		if(p->mgc4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]>1)
        {
			dk = (n+1);
			p->gcorig4[p->mgc4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]-10][5][dk]=1;
        }
	}
    
    
    
//PARA fill
    for(q=0;q<p->gcpara1_count;++q)
    {
    i=p->gcpara1[q][0];
    j=p->gcpara1[q][1];
    k=p->gcpara1[q][2];
        
        if(p->gcpara1[q][6]==1)
        for(n=0;n<p->margin;++n)
        if( p->mgc4[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]>1)
        {
			di = (n+1);
			p->gcorig4[p->mgc4[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]-10][0][di]=1;
        }
    }
    
    for(q=0;q<p->gcpara4_count;++q)
	{
    i=p->gcpara4[q][0];
    j=p->gcpara4[q][1];
    k=p->gcpara4[q][2];
        
        if(p->gcpara4[q][6]==1)
        for(n=0;n<p->margin;++n)
		if( p->mgc4[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]>1)
        {
			di = (n+1);
			p->gcorig4[p->mgc4[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]-10][3][di]=1;
        }
	}

    for(q=0;q<p->gcpara3_count;++q)
    {
    i=p->gcpara3[q][0];
    j=p->gcpara3[q][1];
    k=p->gcpara3[q][2];
        
        if(p->gcpara3[q][6]==1)
        for(n=0;n<p->margin;++n)
        if(p->mgc4[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]>1)
        {
			dj = (n+1);
			p->gcorig4[p->mgc4[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]-10][2][dj]=1;
	    }
    }
    
    for(q=0;q<p->gcpara2_count;++q)
	{
    i=p->gcpara2[q][0];
    j=p->gcpara2[q][1];
    k=p->gcpara2[q][2];
        
        if(p->gcpara2[q][6]==1)
        for(n=0;n<p->margin;++n)
        if( p->mgc4[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]>1)
        {
			dj = (n+1);
			p->gcorig4[p->mgc4[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]-10][1][dj]=1;
        }
	}

	for(q=0;q<p->gcpara5_count;++q)
	{
    i=p->gcpara5[q][0];
    j=p->gcpara5[q][1];
    k=p->gcpara5[q][2];
        
        if(p->gcpara5[q][6]==1)
        for(n=0;n<p->margin;++n)
        if( p->mgc4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]>1)
        {
			dk = (n+1);
			p->gcorig4[p->mgc4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]-10][4][dk]=1;
        }
	}

	for(q=0;q<p->gcpara6_count;++q)
	{
	i=p->gcpara6[q][0];
    j=p->gcpara6[q][1];
    k=p->gcpara6[q][2];
        
        if(p->gcpara6[q][6]==1)
        for(n=0;n<p->margin;++n)
        if(p->mgc4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]>1)
        {
			dk = (n+1);
			p->gcorig4[p->mgc4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]-10][5][dk]=1;
        }
	}
	
}

// -----------------------------------------------------

void mgc4::gcsidefill(lexer *p)
{
	
	fieldint5 side(p);
	int n;
	
	p->Iresize(p->gcside4,p->gcside4_size, p->gcb4_count); 
	p->gcside4_size=p->gcb4_count;
	
	GC4LOOP
	p->gcside4[n]=0;
	
	
	BLOOP
	side(i,j,k)=0;
	
	
	
	GC4LOOP
	{
		i=p->gcb4[n][0];
		j=p->gcb4[n][1];
		k=p->gcb4[n][2];
		
		if(p->gcb4[n][3]==1)
		side(i-1,j,k) +=1;
		
		if(p->gcb4[n][3]==4)
		side(i+1,j,k) +=1;
	}
	
	BLOOP
	if(side(i,j,k)>0)
	side(i,j,k)=10;
	
	GC4LOOP
	{
		i=p->gcb4[n][0];
		j=p->gcb4[n][1];
		k=p->gcb4[n][2];
		
		if(p->gcb4[n][3]==3)
		side(i,j-1,k) +=1;
		
		if(p->gcb4[n][3]==2)
		side(i,j+1,k) +=1;
	}
	
	BLOOP
	{
	if(side(i,j,k)>0 && side(i,j,k)<=2)
	side(i,j,k)=10;
	
	if(side(i,j,k)>10)
	side(i,j,k)=20;
	}
	
	GC4LOOP
	{
		i=p->gcb4[n][0];
		j=p->gcb4[n][1];
		k=p->gcb4[n][2];
		
		if(p->gcb4[n][3]==5)
		side(i,j,k-1) +=1;
		
		if(p->gcb4[n][3]==6)
		side(i,j,k+1) +=1;
	}
	
	BLOOP
	{
	if(side(i,j,k)>0 && side(i,j,k)<=2)
	side(i,j,k)=10;
	
	if(side(i,j,k)>10&& side(i,j,k)<=12)
	side(i,j,k)=20;
	
	if(side(i,j,k)>20)
	side(i,j,k)=30;	
	}
	
	BLOOP
	{
	if(side(i,j,k)==10)
	side(i,j,k)=1;
	
	if(side(i,j,k)==20)
	side(i,j,k)=2;
	
	if(side(i,j,k)==30)
	side(i,j,k)=3;	
	}
	
	
	
	GC4LOOP
	{
		i=p->gcb4[n][0];
		j=p->gcb4[n][1];
		k=p->gcb4[n][2];
		
		if(p->gcb4[n][3]==1)
		p->gcside4[n] = side(i-1,j,k);

		if(p->gcb4[n][3]==4)
		p->gcside4[n] = side(i+1,j,k);
		
		
		if(p->gcb4[n][3]==3)
		p->gcside4[n] = side(i,j-1,k);

		if(p->gcb4[n][3]==2)
		p->gcside4[n] = side(i,j+1,k);
		
		
		if(p->gcb4[n][3]==5)
		p->gcside4[n] = side(i,j,k-1);

		if(p->gcb4[n][3]==6)
		p->gcside4[n] = side(i,j,k+1);
	}
	
}


void mgc4::check_gcb_nbx(lexer *p, ghostcell *gcb)
{
	int gcb4_temp,count;
	
	count=p->gcb4_count;
	gcb4_temp=p->gcb4_count;

	
	
	for(n=0;n<p->gcpara1_count;++n)
    {
    i=p->gcpara1[n][0];
    j=p->gcpara1[n][1];
    k=p->gcpara1[n][2];

        if(p->flag4[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]<0)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]>0)
        ++count;	
    }

    for(n=0;n<p->gcpara2_count;++n)
    {
    i=p->gcpara2[n][0];
    j=p->gcpara2[n][1];
    k=p->gcpara2[n][2];

        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin]<0)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]>0)    
        ++count;		
    }

    for(n=0;n<p->gcpara3_count;++n)
    {
    i=p->gcpara3[n][0];
    j=p->gcpara3[n][1];
    k=p->gcpara3[n][2];

        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]<0)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]>0)
        ++count;
    }

    for(n=0;n<p->gcpara4_count;++n)
    {
    i=p->gcpara4[n][0];
    j=p->gcpara4[n][1];
    k=p->gcpara4[n][2];

        if(p->flag4[(i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]<0)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]>0)
        ++count;
    }

    for(n=0;n<p->gcpara5_count;++n)
    {
    i=p->gcpara5[n][0];
    j=p->gcpara5[n][1];
    k=p->gcpara5[n][2];

        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1]<0)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]>0)
        ++count;
    }

    for(n=0;n<p->gcpara6_count;++n)
    {
    i=p->gcpara6[n][0];
    j=p->gcpara6[n][1];
    k=p->gcpara6[n][2];

        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+1]<0)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]>0)
        ++count;
    }
	
	
	//cout<<p->mpirank<<" COUNT_GCB_GCX   GCB4_old: "<<p->gcb4_count<<" GCB4_new: "<<count<<"   -------- !!!"<<endl;
	
	
	p->Iresize(p->gcb4,p->gcb4_count, count, 6, 6); 
	p->Dresize(p->gcd4,p->gcb4_count, count); 
	
	p->gcb4_count=count;
	
	count=gcb4_temp;
	
		
	
	for(n=0;n<p->gcpara1_count;++n)
    {
    i=p->gcpara1[n][0];
    j=p->gcpara1[n][1];
    k=p->gcpara1[n][2];

        if(p->flag4[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]<0)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]>0)
		{
        p->gcb4[count][0]=i;
		p->gcb4[count][1]=j;
		p->gcb4[count][2]=k;
		p->gcb4[count][3]=1;
		p->gcd4[count]=0.5*p->DXM; 
		
			if(p->flag4[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]==TOPO)
			p->gcb4[count][4]=5;
			
			if(p->flag4[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]==OBJ)
			p->gcb4[count][4]=21;
            
            if(p->flag4[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]==SOLID)
			p->gcb4[count][4]=22;
			
			if(p->flag4[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]==FLT)
			p->gcb4[count][4]=41;
			
		++count;
		}
    }

    for(n=0;n<p->gcpara2_count;++n)
    {
    i=p->gcpara2[n][0];
    j=p->gcpara2[n][1];
    k=p->gcpara2[n][2];

        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin]<0)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]>0)
		{
        p->gcb4[count][0]=i;
		p->gcb4[count][1]=j;
		p->gcb4[count][2]=k;
		p->gcb4[count][3]=2;
		p->gcd4[count]=0.5*p->DXM; 
		
			if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin]==TOPO)
			p->gcb4[count][4]=5;
			
			if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin]==OBJ)
			p->gcb4[count][4]=21;
            
            if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin]==SOLID)
			p->gcb4[count][4]=22;
			
			if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin]==FLT)
			p->gcb4[count][4]=41;
			
		++count;
		}		
    }

    for(n=0;n<p->gcpara3_count;++n)
    {
    i=p->gcpara3[n][0];
    j=p->gcpara3[n][1];
    k=p->gcpara3[n][2];

        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]<0)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]>0)
		{
        p->gcb4[count][0]=i;
		p->gcb4[count][1]=j;
		p->gcb4[count][2]=k;
		p->gcb4[count][3]=3;
		p->gcd4[count]=0.5*p->DXM; 
		
			if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]==TOPO)
			p->gcb4[count][4]=5;
			
			if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]==OBJ)
			p->gcb4[count][4]=21;
            
            if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]==SOLID)
			p->gcb4[count][4]=22;
			
			if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]==FLT)
			p->gcb4[count][4]=41;
			
		++count;
		}
    }

    for(n=0;n<p->gcpara4_count;++n)
    {
    i=p->gcpara4[n][0];
    j=p->gcpara4[n][1];
    k=p->gcpara4[n][2];

        if(p->flag4[(i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]<0)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]>0)
		{
        p->gcb4[count][0]=i;
		p->gcb4[count][1]=j;
		p->gcb4[count][2]=k;
		p->gcb4[count][3]=4;
		p->gcd4[count]=0.5*p->DXM; 
		
			if(p->flag4[(i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]==TOPO)
			p->gcb4[count][4]=5;
			
			if(p->flag4[(i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]==OBJ)
			p->gcb4[count][4]=21;
            
            if(p->flag4[(i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]==SOLID)
			p->gcb4[count][4]=22;
			
			if(p->flag4[(i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]==FLT)
			p->gcb4[count][4]=41;
			
		++count;
		}
    }

    for(n=0;n<p->gcpara5_count;++n)
    {
    i=p->gcpara5[n][0];
    j=p->gcpara5[n][1];
    k=p->gcpara5[n][2];

        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1]<0)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]>0)
		{
        p->gcb4[count][0]=i;
		p->gcb4[count][1]=j;
		p->gcb4[count][2]=k;
		p->gcb4[count][3]=5;
		p->gcd4[count]=0.5*p->DXM; 
		
			if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1]==TOPO)
			p->gcb4[count][4]=5;
			
			if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1]==SOLID)
			p->gcb4[count][4]=21;
            
            if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1]==OBJ)
			p->gcb4[count][4]=22;
			
			if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1]==FLT)
			p->gcb4[count][4]=41;
			
		++count;
		}
    }

    for(n=0;n<p->gcpara6_count;++n)
    {
    i=p->gcpara6[n][0];
    j=p->gcpara6[n][1];
    k=p->gcpara6[n][2];

        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+1]<0)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]>0)
		{
        p->gcb4[count][0]=i;
		p->gcb4[count][1]=j;
		p->gcb4[count][2]=k;
		p->gcb4[count][3]=6;
		p->gcd4[count]=0.5*p->DXM; 
		
			if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+1]==TOPO)
			p->gcb4[count][4]=5;
			
			if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+1]==SOLID)
			p->gcb4[count][4]=21;
            
            if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+1]==OBJ)
			p->gcb4[count][4]=22;
			
			if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+1]==FLT)
			p->gcb4[count][4]=41;
			
		++count;
		}
    }
}



