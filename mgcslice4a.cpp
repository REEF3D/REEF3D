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

#include"mgcslice4a.h"
#include"lexer.h"
#include"ghostcell.h"


mgcslice4a::mgcslice4a(lexer *p)
{
	imin=p->imin;
    imax=p->imax;
    jmin=p->jmin;
    jmax=p->jmax;
}

mgcslice4a::~mgcslice4a()
{
}

void mgcslice4a::makemgc(lexer* p)
{
	p->Iarray(p->mgcsl4a,imax*jmax);

//make gcdir
	p->gcsldirsize4a=1;	
	p->Iarray(p->gcslorig4a, p->gcsldirsize4a, 6,4);
}

void mgcslice4a::mgcsetup(lexer* p)
{
	for(i=0;i<imax*jmax;++i)
	p->mgcsl4a[i]=0;

	SLICELOOP4A
	p->mgcsl4a[IJ]=1;
}

void mgcslice4a::fillmgc(lexer* p)
{
	int q,n;
	
//--------------------------
//WALL1
	QGCSL4ALOOP
	{
        i=p->gcbsl4a[q][0];
        j=p->gcbsl4a[q][1];
		
		if(p->gcbsl4a[q][3]==1)
		for(n=0;n<p->margin;++n)
        p->mgcsl4a[(i-imin-n-1)*jmax + (j-jmin)]+=1;

		if(p->gcbsl4a[q][3]==4)
		for(n=0;n<p->margin;++n)
		p->mgcsl4a[(i-imin+n+1)*jmax + (j-jmin)]+=1;

		if(p->gcbsl4a[q][3]==3)
		for(n=0;n<p->margin;++n)
        p->mgcsl4a[(i-imin)*jmax + (j-jmin-n-1)]+=1;

		if(p->gcbsl4a[q][3]==2)
		for(n=0;n<p->margin;++n)
		p->mgcsl4a[(i-imin)*jmax + (j-jmin+n+1)]+=1;
	}

//--------------------------
//WALL2
    p->gcsl_extra4a=10;
    
	QGCSL4ALOOP
	{
        i=p->gcbsl4a[q][0];
        j=p->gcbsl4a[q][1];

		if(p->gcbsl4a[q][3]==1)
		for(n=0;n<p->margin;++n)
		if(p->mgcsl4a[(i-imin-n-1)*jmax + (j-jmin)]>1
		&& p->mgcsl4a[(i-imin-n-1)*jmax + (j-jmin)]<10)
        {
			p->mgcsl4a[(i-imin-n-1)*jmax + (j-jmin)]=p->gcsl_extra4a;
			++p->gcsl_extra4a;
        }

		if(p->gcbsl4a[q][3]==4)
		for(n=0;n<p->margin;++n)
        if(p->mgcsl4a[(i-imin+n+1)*jmax + (j-jmin)]>1
		&& p->mgcsl4a[(i-imin+n+1)*jmax + (j-jmin)]<10)
        {
			p->mgcsl4a[(i-imin+n+1)*jmax + (j-jmin)]=p->gcsl_extra4a;
			++p->gcsl_extra4a;
        }

		if(p->gcbsl4a[q][3]==3)
		for(n=0;n<p->margin;++n)
		if(p->mgcsl4a[(i-imin)*jmax + (j-jmin-n-1)]>1
		&& p->mgcsl4a[(i-imin)*jmax + (j-jmin-n-1)]<10)
        {
			p->mgcsl4a[(i-imin)*jmax + (j-jmin-n-1)]=p->gcsl_extra4a;
			++p->gcsl_extra4a;
        }

		if(p->gcbsl4a[q][3]==2)
		for(n=0;n<p->margin;++n)
		if(p->mgcsl4a[(i-imin)*jmax + (j-jmin+n+1)]>1
		&& p->mgcsl4a[(i-imin)*jmax + (j-jmin+n+1)]<10)
        {
			p->mgcsl4a[(i-imin)*jmax + (j-jmin+n+1)]=p->gcsl_extra4a;
			++p->gcsl_extra4a;
        }
	}
}

void mgcslice4a::gcdirfill(lexer* p)
{
// GCORIG
    int q,n;
    
	p->Iresize(p->gcslorig4a,p->gcsldirsize4a, p->gcsl_extra4a, 6, 6, 4, 4); 
	p->gcsldirsize4a=p->gcsl_extra4a;
	
	
	for(n=0;n<p->gcsldirsize4a;++n)
	for(q=0;q<6;++q)	
	for(qn=0;qn<4;++qn)	
	p->gcslorig4a[n][q][qn]=0;	
	
	QGCSL4ALOOP
	{
        i=p->gcbsl4a[q][0];
        j=p->gcbsl4a[q][1];

		if(p->gcbsl4a[q][3]==1)
		for(n=0;n<p->margin;++n)
		if( p->mgcsl4a[(i-imin-n-1)*jmax + (j-jmin)]>1)
        {
			di = (n+1);
			p->gcslorig4a[p->mgcsl4a[(i-imin-n-1)*jmax + (j-jmin)]-10][0][di]=1;
        }

		if(p->gcbsl4a[q][3]==4)
		for(n=0;n<p->margin;++n)
        if( p->mgcsl4a[(i-imin+n+1)*jmax + (j-jmin)]>1)
        {
			di = (n+1);
			p->gcslorig4a[p->mgcsl4a[(i-imin+n+1)*jmax + (j-jmin)]-10][3][di]=1;
        }

		if(p->gcbsl4a[q][3]==3)
		for(n=0;n<p->margin;++n)
		if(p->mgcsl4a[(i-imin)*jmax + (j-jmin-n-1)]>1)
        {
			dj = (n+1);
			p->gcslorig4a[p->mgcsl4a[(i-imin)*jmax + (j-jmin-n-1)]-10][2][dj]=1;
	    }

		if(p->gcbsl4a[q][3]==2)
		for(n=0;n<p->margin;++n)
		if( p->mgcsl4a[(i-imin)*jmax + (j-jmin+n+1)]>1)
        {
			dj = (n+1);
			p->gcslorig4a[p->mgcsl4a[(i-imin)*jmax + (j-jmin+n+1)]-10][1][dj]=1;
        }
	}
}

