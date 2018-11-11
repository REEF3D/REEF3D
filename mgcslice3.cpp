/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"mgcslice3.h"
#include"lexer.h"
#include"ghostcell.h"


mgcslice3::mgcslice3(lexer *p)
{
	imin=p->imin;
    imax=p->imax;
    jmin=p->jmin;
    jmax=p->jmax;
}

mgcslice3::~mgcslice3()
{
}

void mgcslice3::makemgc(lexer* p)
{
	p->Iarray(p->mgcsl3,imax*jmax);

//make gcdir
	p->gcsldirsize3=1;	
	p->Iarray(p->gcslorig3, p->gcsldirsize3, 6,4);
}

void mgcslice3::mgcsetup(lexer* p)
{
	for(i=0;i<imax*jmax;++i)
	p->mgcsl3[i]=0;

	SLICELOOP3
	p->mgcsl3[IJ]=1;
}

void mgcslice3::fillmgc(lexer* p)
{
	int q,n;
	
//--------------------------
//WALL1
	QGCSL3LOOP
	{
        i=p->gcbsl3[q][0];
        j=p->gcbsl3[q][1];
		
		if(p->gcbsl3[q][3]==1)
		for(n=0;n<p->margin;++n)
        p->mgcsl3[(i-imin-n-1)*jmax + (j-jmin)]+=1;

		if(p->gcbsl3[q][3]==4)
		for(n=0;n<p->margin;++n)
		p->mgcsl3[(i-imin+n+1)*jmax + (j-jmin)]+=1;

		if(p->gcbsl3[q][3]==3)
		for(n=0;n<p->margin;++n)
        p->mgcsl3[(i-imin)*jmax + (j-jmin-n-1)]+=1;

		if(p->gcbsl3[q][3]==2)
		for(n=0;n<p->margin;++n)
		p->mgcsl3[(i-imin)*jmax + (j-jmin+n+1)]+=1;
	}

//--------------------------
//WALL2
    p->gcsl_extra3=10;
    
	QGCSL3LOOP
	{
        i=p->gcbsl3[q][0];
        j=p->gcbsl3[q][1];

		if(p->gcbsl3[q][3]==1)
		for(n=0;n<p->margin;++n)
		if(p->mgcsl3[(i-imin-n-1)*jmax + (j-jmin)]>1
		&& p->mgcsl3[(i-imin-n-1)*jmax + (j-jmin)]<10)
        {
			p->mgcsl3[(i-imin-n-1)*jmax + (j-jmin)]=p->gcsl_extra3;
			++p->gcsl_extra3;
        }

		if(p->gcbsl3[q][3]==4)
		for(n=0;n<p->margin;++n)
        if(p->mgcsl3[(i-imin+n+1)*jmax + (j-jmin)]>1
		&& p->mgcsl3[(i-imin+n+1)*jmax + (j-jmin)]<10)
        {
			p->mgcsl3[(i-imin+n+1)*jmax + (j-jmin)]=p->gcsl_extra3;
			++p->gcsl_extra3;
        }

		if(p->gcbsl3[q][3]==3)
		for(n=0;n<p->margin;++n)
		if(p->mgcsl3[(i-imin)*jmax + (j-jmin-n-1)]>1
		&& p->mgcsl3[(i-imin)*jmax + (j-jmin-n-1)]<10)
        {
			p->mgcsl3[(i-imin)*jmax + (j-jmin-n-1)]=p->gcsl_extra3;
			++p->gcsl_extra3;
        }

		if(p->gcbsl3[q][3]==2)
		for(n=0;n<p->margin;++n)
		if(p->mgcsl3[(i-imin)*jmax + (j-jmin+n+1)]>1
		&& p->mgcsl3[(i-imin)*jmax + (j-jmin+n+1)]<10)
        {
			p->mgcsl3[(i-imin)*jmax + (j-jmin+n+1)]=p->gcsl_extra3;
			++p->gcsl_extra3;
        }
	}
}

void mgcslice3::gcdirfill(lexer* p)
{
// GCORIG
    int q,n;
    
	p->Iresize(p->gcslorig3,p->gcsldirsize3, p->gcsl_extra3, 6, 6, 4, 4); 
	p->gcsldirsize3=p->gcsl_extra3;
	
	
	for(n=0;n<p->gcsldirsize3;++n)
	for(q=0;q<6;++q)	
	for(qn=0;qn<4;++qn)	
	p->gcslorig3[n][q][qn]=0;	
	
	QGCSL3LOOP
	{
        i=p->gcbsl3[q][0];
        j=p->gcbsl3[q][1];

		if(p->gcbsl3[q][3]==1)
		for(n=0;n<p->margin;++n)
		if( p->mgcsl3[(i-imin-n-1)*jmax + (j-jmin)]>1)
        {
			di = (n+1);
			p->gcslorig3[p->mgcsl3[(i-imin-n-1)*jmax + (j-jmin)]-10][0][di]=1;
        }

		if(p->gcbsl3[q][3]==4)
		for(n=0;n<p->margin;++n)
        if( p->mgcsl3[(i-imin+n+1)*jmax + (j-jmin)]>1)
        {
			di = (n+1);
			p->gcslorig3[p->mgcsl3[(i-imin+n+1)*jmax + (j-jmin)]-10][3][di]=1;
        }

		if(p->gcbsl3[q][3]==3)
		for(n=0;n<p->margin;++n)
		if(p->mgcsl3[(i-imin)*jmax + (j-jmin-n-1)]>1)
        {
			dj = (n+1);
			p->gcslorig3[p->mgcsl3[(i-imin)*jmax + (j-jmin-n-1)]-10][2][dj]=1;
	    }

		if(p->gcbsl3[q][3]==2)
		for(n=0;n<p->margin;++n)
		if( p->mgcsl3[(i-imin)*jmax + (j-jmin+n+1)]>1)
        {
			dj = (n+1);
			p->gcslorig3[p->mgcsl3[(i-imin)*jmax + (j-jmin+n+1)]-10][1][dj]=1;
        }
	}
}

