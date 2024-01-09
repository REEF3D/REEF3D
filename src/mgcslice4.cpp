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

#include"mgcslice4.h"
#include"lexer.h"
#include"ghostcell.h"


mgcslice4::mgcslice4(lexer *p)
{
	imin=p->imin;
    imax=p->imax;
    jmin=p->jmin;
    jmax=p->jmax;
}

mgcslice4::~mgcslice4()
{
}

void mgcslice4::makemgc(lexer* p)
{
	p->Iarray(p->mgcsl4,imax*jmax);

//make gcdir
	p->gcsldirsize4=1;	
	p->Iarray(p->gcslorig4, p->gcsldirsize4, 6,4);

}

void mgcslice4::mgcsetup(lexer* p)
{
	for(i=0;i<imax*jmax;++i)
	p->mgcsl4[i]=0;

	SLICELOOP4
	p->mgcsl4[IJ]=1;
}

void mgcslice4::fillmgc(lexer* p)
{
	int q,n;
	
//--------------------------
//WALL1
	QGCSL4LOOP
	{
        i=p->gcbsl4[q][0];
        j=p->gcbsl4[q][1];
		
		if(p->gcbsl4[q][3]==1)
		for(n=0;n<p->margin;++n)
        p->mgcsl4[(i-imin-n-1)*jmax + (j-jmin)]+=1;

		if(p->gcbsl4[q][3]==4)
		for(n=0;n<p->margin;++n)
		p->mgcsl4[(i-imin+n+1)*jmax + (j-jmin)]+=1;

		if(p->gcbsl4[q][3]==3)
		for(n=0;n<p->margin;++n)
        p->mgcsl4[(i-imin)*jmax + (j-jmin-n-1)]+=1;

		if(p->gcbsl4[q][3]==2)
		for(n=0;n<p->margin;++n)
		p->mgcsl4[(i-imin)*jmax + (j-jmin+n+1)]+=1;
	}

//--------------------------
//WALL2
    p->gcsl_extra4=10;
    
	QGCSL4LOOP
	{
        i=p->gcbsl4[q][0];
        j=p->gcbsl4[q][1];

		if(p->gcbsl4[q][3]==1)
		for(n=0;n<p->margin;++n)
		if(p->mgcsl4[(i-imin-n-1)*jmax + (j-jmin)]>1
		&& p->mgcsl4[(i-imin-n-1)*jmax + (j-jmin)]<10)
        {
			p->mgcsl4[(i-imin-n-1)*jmax + (j-jmin)]=p->gcsl_extra4;
			++p->gcsl_extra4;
        }

		if(p->gcbsl4[q][3]==4)
		for(n=0;n<p->margin;++n)
        if(p->mgcsl4[(i-imin+n+1)*jmax + (j-jmin)]>1
		&& p->mgcsl4[(i-imin+n+1)*jmax + (j-jmin)]<10)
        {
			p->mgcsl4[(i-imin+n+1)*jmax + (j-jmin)]=p->gcsl_extra4;
			++p->gcsl_extra4;
        }

		if(p->gcbsl4[q][3]==3)
		for(n=0;n<p->margin;++n)
		if(p->mgcsl4[(i-imin)*jmax + (j-jmin-n-1)]>1
		&& p->mgcsl4[(i-imin)*jmax + (j-jmin-n-1)]<10)
        {
			p->mgcsl4[(i-imin)*jmax + (j-jmin-n-1)]=p->gcsl_extra4;
			++p->gcsl_extra4;
        }

		if(p->gcbsl4[q][3]==2)
		for(n=0;n<p->margin;++n)
		if(p->mgcsl4[(i-imin)*jmax + (j-jmin+n+1)]>1
		&& p->mgcsl4[(i-imin)*jmax + (j-jmin+n+1)]<10)
        {
			p->mgcsl4[(i-imin)*jmax + (j-jmin+n+1)]=p->gcsl_extra4;
			++p->gcsl_extra4;
        }
	}	    
}

void mgcslice4::gcdirfill(lexer* p)
{
// GCORIG
    int q,n;
    
	p->Iresize(p->gcslorig4,p->gcsldirsize4, p->gcsl_extra4, 6, 6, 4, 4); 
	p->gcsldirsize4=p->gcsl_extra4;
	
	
	for(n=0;n<p->gcsldirsize4;++n)
	for(q=0;q<6;++q)	
	for(qn=0;qn<4;++qn)	
	p->gcslorig4[n][q][qn]=0;	
	
	QGCSL4LOOP
	{
        i=p->gcbsl4[q][0];
        j=p->gcbsl4[q][1];

		if(p->gcbsl4[q][3]==1)
		for(n=0;n<p->margin;++n)
		if( p->mgcsl4[(i-imin-n-1)*jmax + (j-jmin)]>1)
        {
			di = (n+1);
			p->gcslorig4[p->mgcsl4[(i-imin-n-1)*jmax + (j-jmin)]-10][0][di]=1;
        }

		if(p->gcbsl4[q][3]==4)
		for(n=0;n<p->margin;++n)
        if( p->mgcsl4[(i-imin+n+1)*jmax + (j-jmin)]>1)
        {
			di = (n+1);
			p->gcslorig4[p->mgcsl4[(i-imin+n+1)*jmax + (j-jmin)]-10][3][di]=1;
        }

		if(p->gcbsl4[q][3]==3)
		for(n=0;n<p->margin;++n)
		if(p->mgcsl4[(i-imin)*jmax + (j-jmin-n-1)]>1)
        {
			dj = (n+1);
			p->gcslorig4[p->mgcsl4[(i-imin)*jmax + (j-jmin-n-1)]-10][2][dj]=1;
	    }

		if(p->gcbsl4[q][3]==2)
		for(n=0;n<p->margin;++n)
		if( p->mgcsl4[(i-imin)*jmax + (j-jmin+n+1)]>1)
        {
			dj = (n+1);
			p->gcslorig4[p->mgcsl4[(i-imin)*jmax + (j-jmin+n+1)]-10][1][dj]=1;
        }
	}
}

