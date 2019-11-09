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

#include"mgcslice3.h"
#include"cart3.h"
#include"lexer.h"
#include"fdm.h"

void mgcslice3::make_ggc(lexer* p)
{
	p->ggcslsize3=1;
	p->Iarray(p->ggcsl3,p->ggcslsize3,3);
}

void mgcslice3::fill_ggc(lexer* p)
{
	int q,qq,n,nn,a;
	int check;
	
	p->Iarray(p->ggcslmem3,imax*jmax);

//--------------------------
//WALL1

	GCSL3LOOP
	{
        i=p->gcbsl3[n][0];
		j=p->gcbsl3[n][1];

		if(p->gcbsl3[n][3]==1)
		for(q=0;q<p->margin;++q)
        p->ggcslmem3[(i-imin-q-1)*jmax + (j-jmin)]+=1;

		if(p->gcbsl3[n][3]==4)
		for(q=0;q<p->margin;++q)
		p->ggcslmem3[(i-imin+q+1)*jmax + (j-jmin)]+=1;

		if(p->gcbsl3[n][3]==3)
		for(q=0;q<p->margin;++q)
        p->ggcslmem3[(i-imin)*jmax + (j-jmin-q-1)]+=1;

		if(p->gcbsl3[n][3]==2)
		for(q=0;q<p->margin;++q)
		p->ggcslmem3[(i-imin)*jmax + (j-jmin+q+1)]+=1;

	}

// count entries
	p->ggcslcount3=0;
	a=0;
	for(i=0;i<imax;++i)
	for(j=0;j<jmax;++j)
	{
        if(p->ggcslmem3[a]>1)
        ++p->ggcslcount3;

	++a;
	}
	
	p->Iresize(p->ggcsl3,p->ggcslsize3,p->ggcslcount3*p->margin, 3, 3);
	p->ggcslsize3=p->ggcslcount3*p->margin;

//--------------------------
//WALL2
	n=0;
	QQGCSL3LOOP
	{
        i=p->gcbsl3[qq][0];
		j=p->gcbsl3[qq][1];

		if(p->gcbsl3[qq][3]==1)
		for(q=0;q<p->margin;++q)
		if(p->ggcslmem3[(i-imin-q-1)*jmax + (j-jmin)]>1)
        {
            if(p->ggcslmem3[(i-imin-q-1)*jmax + (j-jmin)]<10)
            {
             p->ggcslmem3[(i-imin-q-1)*jmax + (j-jmin)]=n+10;
			 p->ggcsl3[n][0]=i-q-1;
			 p->ggcsl3[n][1]=j;
			 ++n;
            }
        }

        if(p->gcbsl3[qq][3]==4)
        for(q=0;q<p->margin;++q)
		if(p->ggcslmem3[(i-imin+q+1)*jmax + (j-jmin)]>1)
        {
            if(p->ggcslmem3[(i-imin+q+1)*jmax + (j-jmin)]<10)
            {
            p->ggcslmem3[(i-imin+q+1)*jmax + (j-jmin)]=n+10;
			p->ggcsl3[n][0]=i+q+1;
			p->ggcsl3[n][1]=j;
			++n;
            }
        }

        if(p->gcbsl3[qq][3]==3)
        for(q=0;q<p->margin;++q)
		if(p->ggcslmem3[(i-imin)*jmax + (j-jmin-q-1)]>1)
        {
            if(p->ggcslmem3[(i-imin)*jmax + (j-jmin-q-1)]<10)
            {
            p->ggcslmem3[(i-imin)*jmax + (j-jmin-q-1)]=n+10;
			p->ggcsl3[n][0]=i;
			p->ggcsl3[n][1]=j-q-1;
			++n;
            }
        }

		if(p->gcbsl3[qq][3]==2)
		for(q=0;q<p->margin;++q)
		if(p->ggcslmem3[(i-imin)*jmax + (j-jmin+q+1)]>1)
        {
            if(p->ggcslmem3[(i-imin)*jmax + (j-jmin+q+1)]<10)
            {
            p->ggcslmem3[(i-imin)*jmax + (j-jmin+q+1)]=n+10;
			p->ggcsl3[n][0]=i;
			p->ggcsl3[n][1]=j+q+1;
			++n;
            }
        }
	}
	p->ggcslcount3=n;


	p->del_Iarray(p->ggcslmem3,imax*jmax);
}


