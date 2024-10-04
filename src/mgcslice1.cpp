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

#include"mgcslice1.h"
#include"lexer.h"
#include"ghostcell.h"


mgcslice1::mgcslice1(lexer *p)
{
	imin=p->imin;
    imax=p->imax;
    jmin=p->jmin;
    jmax=p->jmax;
}

mgcslice1::~mgcslice1()
{
}

void mgcslice1::makemgc(lexer* p)
{
	p->Iarray(p->mgcsl1,imax*jmax);

//make gcdir
	p->gcsldirsize1=1;	
	p->Iarray(p->gcslorig1, p->gcsldirsize1, 6,4);
    
//flag1
    for(i=0;i<p->imax*p->jmax; ++i)
	{
	p->flagslice1[i]=p->flagslice4[i];
	}

    SLICELOOP4
    {

        if(p->flagslice4[Ip1J]<0)
        p->flagslice1[IJ]=p->flagslice4[Ip1J];
    }
}

void mgcslice1::mgcsetup(lexer* p)
{
	for(i=0;i<imax*jmax;++i)
	p->mgcsl1[i]=0;

	SLICELOOP1
	p->mgcsl1[IJ]=1;
}

void mgcslice1::fillmgc(lexer* p)
{
	int q,n;
	
//--------------------------
//WALL1
	QGCSL1LOOP
	{
        i=p->gcbsl1[q][0];
        j=p->gcbsl1[q][1];
		
		if(p->gcbsl1[q][3]==1)
		for(n=0;n<p->margin;++n)
        p->mgcsl1[(i-imin-n-1)*jmax + (j-jmin)]+=1;

		if(p->gcbsl1[q][3]==4)
		for(n=0;n<p->margin;++n)
		p->mgcsl1[(i-imin+n+1)*jmax + (j-jmin)]+=1;

		if(p->gcbsl1[q][3]==3)
		for(n=0;n<p->margin;++n)
        p->mgcsl1[(i-imin)*jmax + (j-jmin-n-1)]+=1;

		if(p->gcbsl1[q][3]==2)
		for(n=0;n<p->margin;++n)
		p->mgcsl1[(i-imin)*jmax + (j-jmin+n+1)]+=1;
	}

//--------------------------
//WALL2
    p->gcsl_extra1=10;
    
	QGCSL1LOOP
	{
        i=p->gcbsl1[q][0];
        j=p->gcbsl1[q][1];

		if(p->gcbsl1[q][3]==1)
		for(n=0;n<p->margin;++n)
		if(p->mgcsl1[(i-imin-n-1)*jmax + (j-jmin)]>1
		&& p->mgcsl1[(i-imin-n-1)*jmax + (j-jmin)]<10)
        {
			p->mgcsl1[(i-imin-n-1)*jmax + (j-jmin)]=p->gcsl_extra1;
			++p->gcsl_extra1;
        }

		if(p->gcbsl1[q][3]==4)
		for(n=0;n<p->margin;++n)
        if(p->mgcsl1[(i-imin+n+1)*jmax + (j-jmin)]>1
		&& p->mgcsl1[(i-imin+n+1)*jmax + (j-jmin)]<10)
        {
			p->mgcsl1[(i-imin+n+1)*jmax + (j-jmin)]=p->gcsl_extra1;
			++p->gcsl_extra1;
        }

		if(p->gcbsl1[q][3]==3)
		for(n=0;n<p->margin;++n)
		if(p->mgcsl1[(i-imin)*jmax + (j-jmin-n-1)]>1
		&& p->mgcsl1[(i-imin)*jmax + (j-jmin-n-1)]<10)
        {
			p->mgcsl1[(i-imin)*jmax + (j-jmin-n-1)]=p->gcsl_extra1;
			++p->gcsl_extra1;
        }

		if(p->gcbsl1[q][3]==2)
		for(n=0;n<p->margin;++n)
		if(p->mgcsl1[(i-imin)*jmax + (j-jmin+n+1)]>1
		&& p->mgcsl1[(i-imin)*jmax + (j-jmin+n+1)]<10)
        {
			p->mgcsl1[(i-imin)*jmax + (j-jmin+n+1)]=p->gcsl_extra1;
			++p->gcsl_extra1;
        }
	}
}

void mgcslice1::gcdirfill(lexer* p)
{
// GCORIG
    int q,n;
    
	p->Iresize(p->gcslorig1,p->gcsldirsize1, p->gcsl_extra1, 6, 6, 4, 4); 
	p->gcsldirsize1=p->gcsl_extra1;
	
	
	for(n=0;n<p->gcsldirsize1;++n)
	for(q=0;q<6;++q)	
	for(qn=0;qn<4;++qn)	
	p->gcslorig1[n][q][qn]=0;	
	
	QGCSL1LOOP
	{
        i=p->gcbsl1[q][0];
        j=p->gcbsl1[q][1];

		if(p->gcbsl1[q][3]==1)
		for(n=0;n<p->margin;++n)
		if( p->mgcsl1[(i-imin-n-1)*jmax + (j-jmin)]>1)
        {
			di = (n+1);
			p->gcslorig1[p->mgcsl1[(i-imin-n-1)*jmax + (j-jmin)]-10][0][di]=1;
        }

		if(p->gcbsl1[q][3]==4)
		for(n=0;n<p->margin;++n)
        if( p->mgcsl1[(i-imin+n+1)*jmax + (j-jmin)]>1)
        {
			di = (n+1);
			p->gcslorig1[p->mgcsl1[(i-imin+n+1)*jmax + (j-jmin)]-10][3][di]=1;
        }

		if(p->gcbsl1[q][3]==3)
		for(n=0;n<p->margin;++n)
		if(p->mgcsl1[(i-imin)*jmax + (j-jmin-n-1)]>1)
        {
			dj = (n+1);
			p->gcslorig1[p->mgcsl1[(i-imin)*jmax + (j-jmin-n-1)]-10][2][dj]=1;
	    }

		if(p->gcbsl1[q][3]==2)
		for(n=0;n<p->margin;++n)
		if( p->mgcsl1[(i-imin)*jmax + (j-jmin+n+1)]>1)
        {
			dj = (n+1);
			p->gcslorig1[p->mgcsl1[(i-imin)*jmax + (j-jmin+n+1)]-10][1][dj]=1;
        }
	}
}

