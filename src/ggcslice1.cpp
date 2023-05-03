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
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"mgcslice1.h"
#include"cart1.h"
#include"lexer.h"
#include"fdm.h"

void mgcslice1::make_ggc(lexer* p)
{
	p->ggcslsize1=1;
	p->Iarray(p->ggcsl1,p->ggcslsize1,3);
}

void mgcslice1::fill_ggc(lexer* p)
{
	int q,qq,n,nn,a;
	int check;

	p->Iarray(p->ggcslmem1,imax*jmax);

//--------------------------
//WALL1

	GCSL1LOOP
	{
        i=p->gcbsl1[n][0];
		j=p->gcbsl1[n][1];

		if(p->gcbsl1[n][3]==1)
		for(q=0;q<p->margin;++q)
        p->ggcslmem1[(i-imin-q-1)*jmax + (j-jmin)]+=1;

		if(p->gcbsl1[n][3]==4)
		for(q=0;q<p->margin;++q)
		p->ggcslmem1[(i-imin+q+1)*jmax + (j-jmin)]+=1;

		if(p->gcbsl1[n][3]==3)
		for(q=0;q<p->margin;++q)
        p->ggcslmem1[(i-imin)*jmax + (j-jmin-q-1)]+=1;

		if(p->gcbsl1[n][3]==2)
		for(q=0;q<p->margin;++q)
		p->ggcslmem1[(i-imin)*jmax + (j-jmin+q+1)]+=1;

	}

// count entries
	p->ggcslcount1=0;
	a=0;
	for(i=0;i<imax;++i)
	for(j=0;j<jmax;++j)
	{
        if(p->ggcslmem1[a]>1)
        ++p->ggcslcount1;

	++a;
	}

	p->Iresize(p->ggcsl1,p->ggcslsize1,p->ggcslcount1*p->margin, 3, 3);
	p->ggcslsize1=p->ggcslcount1*p->margin;

//--------------------------
//WALL2
	n=0;
	QQGCSL1LOOP
	{
        i=p->gcbsl1[qq][0];
		j=p->gcbsl1[qq][1];

		if(p->gcbsl1[qq][3]==1)
		for(q=0;q<p->margin;++q)
		if(p->ggcslmem1[(i-imin-q-1)*jmax + (j-jmin)]>1)
        {
            if(p->ggcslmem1[(i-imin-q-1)*jmax + (j-jmin)]<10)
            {
             p->ggcslmem1[(i-imin-q-1)*jmax + (j-jmin)]=n+10;
			 p->ggcsl1[n][0]=i-q-1;
			 p->ggcsl1[n][1]=j;
			 ++n;
            }
        }

        if(p->gcbsl1[qq][3]==4)
        for(q=0;q<p->margin;++q)
		if(p->ggcslmem1[(i-imin+q+1)*jmax + (j-jmin)]>1)
        {
            if(p->ggcslmem1[(i-imin+q+1)*jmax + (j-jmin)]<10)
            {
            p->ggcslmem1[(i-imin+q+1)*jmax + (j-jmin)]=n+10;
			p->ggcsl1[n][0]=i+q+1;
			p->ggcsl1[n][1]=j;
			++n;
            }
        }

        if(p->gcbsl1[qq][3]==3)
        for(q=0;q<p->margin;++q)
		if(p->ggcslmem1[(i-imin)*jmax + (j-jmin-q-1)]>1)
        {
            if(p->ggcslmem1[(i-imin)*jmax + (j-jmin-q-1)]<10)
            {
            p->ggcslmem1[(i-imin)*jmax + (j-jmin-q-1)]=n+10;
			p->ggcsl1[n][0]=i;
			p->ggcsl1[n][1]=j-q-1;
			++n;
            }
        }

		if(p->gcbsl1[qq][3]==2)
		for(q=0;q<p->margin;++q)
		if(p->ggcslmem1[(i-imin)*jmax + (j-jmin+q+1)]>1)
        {
            if(p->ggcslmem1[(i-imin)*jmax + (j-jmin+q+1)]<10)
            {
            p->ggcslmem1[(i-imin)*jmax + (j-jmin+q+1)]=n+10;
			p->ggcsl1[n][0]=i;
			p->ggcsl1[n][1]=j+q+1;
			++n;
            }
        }
	}
	p->ggcslcount1=n;

	p->del_Iarray(p->ggcslmem1,imax*jmax);
}


