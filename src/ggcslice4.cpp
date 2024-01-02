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
#include"cart4.h"
#include"lexer.h"
#include"fdm.h"

void mgcslice4::make_ggc(lexer* p)
{
	p->ggcslsize4=1;
	p->Iarray(p->ggcsl4,p->ggcslsize4,3);
}

void mgcslice4::fill_ggc(lexer* p)
{
	int q,qq,n,nn,a;
	int check;
	
	p->Iarray(p->ggcslmem4,imax*jmax);

//--------------------------
//WALL1

	GCSL4LOOP
	{
        i=p->gcbsl4[n][0];
		j=p->gcbsl4[n][1];

		if(p->gcbsl4[n][3]==1)
		for(q=0;q<p->margin;++q)
        p->ggcslmem4[(i-imin-q-1)*jmax + (j-jmin)]+=1;

		if(p->gcbsl4[n][3]==4)
		for(q=0;q<p->margin;++q)
		p->ggcslmem4[(i-imin+q+1)*jmax + (j-jmin)]+=1;

		if(p->gcbsl4[n][3]==3)
		for(q=0;q<p->margin;++q)
        p->ggcslmem4[(i-imin)*jmax + (j-jmin-q-1)]+=1;

		if(p->gcbsl4[n][3]==2)
		for(q=0;q<p->margin;++q)
		p->ggcslmem4[(i-imin)*jmax + (j-jmin+q+1)]+=1;

	}

// count entries
	p->ggcslcount4=0;
	a=0;
	for(i=0;i<imax;++i)
	for(j=0;j<jmax;++j)
	{
        if(p->ggcslmem4[a]>1)
        ++p->ggcslcount4;

	++a;
	}
	
	p->Iresize(p->ggcsl4,p->ggcslsize4,p->ggcslcount4*p->margin, 3, 3);
	p->ggcslsize4=p->ggcslcount4*p->margin;

//--------------------------
//WALL2
	n=0;
	QQGCSL4LOOP
	{
        i=p->gcbsl4[qq][0];
		j=p->gcbsl4[qq][1];

		if(p->gcbsl4[qq][3]==1)
		for(q=0;q<p->margin;++q)
		if(p->ggcslmem4[(i-imin-q-1)*jmax + (j-jmin)]>1)
        {
            if(p->ggcslmem4[(i-imin-q-1)*jmax + (j-jmin)]<10)
            {
             p->ggcslmem4[(i-imin-q-1)*jmax + (j-jmin)]=n+10;
			 p->ggcsl4[n][0]=i-q-1;
			 p->ggcsl4[n][1]=j;
			 ++n;
            }
        }

        if(p->gcbsl4[qq][3]==4)
        for(q=0;q<p->margin;++q)
		if(p->ggcslmem4[(i-imin+q+1)*jmax + (j-jmin)]>1)
        {
            if(p->ggcslmem4[(i-imin+q+1)*jmax + (j-jmin)]<10)
            {
            p->ggcslmem4[(i-imin+q+1)*jmax + (j-jmin)]=n+10;
			p->ggcsl4[n][0]=i+q+1;
			p->ggcsl4[n][1]=j;
			++n;
            }
        }

        if(p->gcbsl4[qq][3]==3)
        for(q=0;q<p->margin;++q)
		if(p->ggcslmem4[(i-imin)*jmax + (j-jmin-q-1)]>1)
        {
            if(p->ggcslmem4[(i-imin)*jmax + (j-jmin-q-1)]<10)
            {
            p->ggcslmem4[(i-imin)*jmax + (j-jmin-q-1)]=n+10;
			p->ggcsl4[n][0]=i;
			p->ggcsl4[n][1]=j-q-1;
			++n;
            }
        }

		if(p->gcbsl4[qq][3]==2)
		for(q=0;q<p->margin;++q)
		if(p->ggcslmem4[(i-imin)*jmax + (j-jmin+q+1)]>1)
        {
            if(p->ggcslmem4[(i-imin)*jmax + (j-jmin+q+1)]<10)
            {
            p->ggcslmem4[(i-imin)*jmax + (j-jmin+q+1)]=n+10;
			p->ggcsl4[n][0]=i;
			p->ggcsl4[n][1]=j+q+1;
			++n;
            }
        }
	}
	p->ggcslcount4=n;

    
    //cout<<p->mpirank<<"  ggcslcount: "<<p->ggcslcount4<<endl;
	p->del_Iarray(p->ggcslmem4,imax*jmax);
}


