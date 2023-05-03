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

#include"mgc6.h"
#include"cart4.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fieldint5.h"

mgc6::mgc6(lexer *p)
{
	imin=p->imin;
    imax=p->imax;
    jmin=p->jmin;
    jmax=p->jmax;
    kmin=p->kmin;
    kmax=p->kmax;
}

mgc6::~mgc6()
{
}

void mgc6::makemgc(lexer* p)
{
	p->Iarray(p->mgc6,kmax*jmax*imax);

//make gcdir
	p->gcdirsize6=1;	
	p->Iarray(p->gcorig6, p->gcdirsize6, 6,4);
}

void mgc6::mgcsetup(lexer* p)
{
	for(i=0;i<imax*jmax*kmax;++i)
	p->mgc6[i]=0;

	BASELOOP
	p->mgc6[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin]=1;
}

void mgc6::fillmgc(lexer* p)
{
	int q,n;
	
//--------------------------
//WALL1
	QGC6LOOP
	{
        i=p->gcb4[q][0];
        j=p->gcb4[q][1];
        k=p->gcb4[q][2];
		
		if(p->gcb4[q][3]==1)
		for(n=0;n<p->margin;++n)
        p->mgc6[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]+=1;

		if(p->gcb4[q][3]==4)
		for(n=0;n<p->margin;++n)
		p->mgc6[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]+=1;

		if(p->gcb4[q][3]==3)
		for(n=0;n<p->margin;++n)
        p->mgc6[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]+=1;

		if(p->gcb4[q][3]==2)
		for(n=0;n<p->margin;++n)
		p->mgc6[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]+=1;

		if(p->gcb4[q][3]==5)
		for(n=0;n<p->margin;++n)
		p->mgc6[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]+=1;

		if(p->gcb4[q][3]==6)
		for(n=0;n<p->margin;++n)
		p->mgc6[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]+=1;
	}
    
//PARA1
    for(q=0;q<p->gcpara1_count;++q)
    {
    i=p->gcpara1[q][0];
    j=p->gcpara1[q][1];
    k=p->gcpara1[q][2];
        
        if(p->gcpara1[q][8]==1)
        for(n=0;n<p->margin;++n)
        p->mgc6[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]+=1;
    }
    
    for(q=0;q<p->gcpara4_count;++q)
	{
    i=p->gcpara4[q][0];
    j=p->gcpara4[q][1];
    k=p->gcpara4[q][2];
        
        if(p->gcpara4[q][8]==1)
        for(n=0;n<p->margin;++n)
		p->mgc6[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]+=1;
	}

    for(q=0;q<p->gcpara3_count;++q)
    {
    i=p->gcpara3[q][0];
    j=p->gcpara3[q][1];
    k=p->gcpara3[q][2];
        
        if(p->gcpara3[q][8]==1)
        for(n=0;n<p->margin;++n)
        p->mgc6[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]+=1;
    }
    
    for(q=0;q<p->gcpara2_count;++q)
	{
    i=p->gcpara2[q][0];
    j=p->gcpara2[q][1];
    k=p->gcpara2[q][2];
        
        if(p->gcpara2[q][8]==1)
        for(n=0;n<p->margin;++n)
        p->mgc6[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]+=1;
	}

	for(q=0;q<p->gcpara5_count;++q)
	{
    i=p->gcpara5[q][0];
    j=p->gcpara5[q][1];
    k=p->gcpara5[q][2];
        
        if(p->gcpara5[q][8]==1)
        for(n=0;n<p->margin;++n)
        p->mgc6[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]+=1;
	}

	for(q=0;q<p->gcpara6_count;++q)
	{
	i=p->gcpara6[q][0];
    j=p->gcpara6[q][1];
    k=p->gcpara6[q][2];
        
        if(p->gcpara6[q][8]==1)
        for(n=0;n<p->margin;++n)
        p->mgc6[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]+=1;
	}

//--------------------------
//WALL2
    p->gcextra6=10;
    
	QGC6LOOP
	{
        i=p->gcb4[q][0];
        j=p->gcb4[q][1];
        k=p->gcb4[q][2];

		if(p->gcb4[q][3]==1)
		for(n=0;n<p->margin;++n)
		if(p->mgc6[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]>1
		&& p->mgc6[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]<10)
        {
			p->mgc6[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]=p->gcextra6;
			++p->gcextra6;
        }

		if(p->gcb4[q][3]==4)
		for(n=0;n<p->margin;++n)
        if(p->mgc6[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]>1
		&& p->mgc6[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]<10)
        {
			p->mgc6[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]=p->gcextra6;
			++p->gcextra6;
        }

		if(p->gcb4[q][3]==3)
		for(n=0;n<p->margin;++n)
		if(p->mgc6[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]>1
		&& p->mgc6[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]<10)
        {
			p->mgc6[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]=p->gcextra6;
			++p->gcextra6;
        }

		if(p->gcb4[q][3]==2)
		for(n=0;n<p->margin;++n)
		if(p->mgc6[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]>1
		&& p->mgc6[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]<10)
        {
			p->mgc6[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]=p->gcextra6;
			++p->gcextra6;
        }

		if(p->gcb4[q][3]==5)
		for(n=0;n<p->margin;++n)
		if( p->mgc6[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]>1
		&& p->mgc6[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]<10)
        {
			p->mgc6[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]=p->gcextra6;
			++p->gcextra6;
        }

		if(p->gcb4[q][3]==6)
		for(n=0;n<p->margin;++n)
		if(p->mgc6[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]>1
		&& p->mgc6[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]<10)
        {
			p->mgc6[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]=p->gcextra6;
			++p->gcextra6;
        }
	}
    
}

void mgc6::gcdirfill(lexer* p)
{
// GCORIG
    int q,n;
    
	p->Iresize(p->gcorig6,p->gcdirsize6, p->gcextra6, 6, 6, 4, 4); 
	p->gcdirsize6=p->gcextra6;
	
	
	for(n=0;n<p->gcdirsize6;++n)
	for(q=0;q<6;++q)	
	for(qn=0;qn<4;++qn)	
	p->gcorig6[n][q][qn]=0;	
	
	QGC6LOOP
	{
        i=p->gcb4[q][0];
        j=p->gcb4[q][1];
        k=p->gcb4[q][2];

		if(p->gcb4[q][3]==1)
		for(n=0;n<p->margin;++n)
		if( p->mgc6[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]>1)
        {
			di = (n+1);
			p->gcorig6[p->mgc6[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]-10][0][di]=1;
        }

		if(p->gcb4[q][3]==4)
		for(n=0;n<p->margin;++n)
        if( p->mgc6[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]>1)
        {
			di = (n+1);
			p->gcorig6[p->mgc6[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]-10][3][di]=1;
        }

		if(p->gcb4[q][3]==3)
		for(n=0;n<p->margin;++n)
		if(p->mgc6[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]>1)
        {
			dj = (n+1);
			p->gcorig6[p->mgc6[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]-10][2][dj]=1;
	    }

		if(p->gcb4[q][3]==2)
		for(n=0;n<p->margin;++n)
		if( p->mgc6[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]>1)
        {
			dj = (n+1);
			p->gcorig6[p->mgc6[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]-10][1][dj]=1;
        }

		if(p->gcb4[q][3]==5)
		for(n=0;n<p->margin;++n)
		if( p->mgc6[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]>1)
        {
			dk = (n+1);
			p->gcorig6[p->mgc6[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]-10][4][dk]=1;
        }

		if(p->gcb4[q][3]==6)
		for(n=0;n<p->margin;++n)
		if(p->mgc6[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]>1)
        {
			dk = (n+1);
			p->gcorig6[p->mgc6[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]-10][5][dk]=1;
        }
	}
    
}
