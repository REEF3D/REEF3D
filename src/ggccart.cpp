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

#include"ggccart.h"
#include"lexer.h"
#include"fdm.h"

void ggccart::make_ggc(lexer* p, int ggcsize, int **ggc)
{
	ggcsize=1;
	p->Iarray(ggc,ggcsize,12);
}

void ggccart::startggc(lexer* p)
{
}

void ggccart::fill_ggc(lexer* p, int gcb_count, int **gcb, int ggccount, int **ggc)
{
	int q,qq,n,nn,a;
	int check;
	
	p->Iarray(ggcmem,kmax*jmax*imax);

//--------------------------
//WALL1

	GC4LOOP
	{
       i=p->gcb4[n][0];
		j=p->gcb4[n][1];
		k=p->gcb4[n][2];

		if(p->gcb4[n][3]==1)
		for(q=0;q<p->margin;++q)
		if(p->flag4[(i-imin-q-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]>0)
        ggcmem[(i-imin-q-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]+=1;

		if(p->gcb4[n][3]==4)
		for(q=0;q<p->margin;++q)
		if(p->flag4[(i-imin+q+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]>0)
		ggcmem[(i-imin+q+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]+=1;

		if(p->gcb4[n][3]==3)
		for(q=0;q<p->margin;++q)
		if(p->flag4[(i-imin)*jmax*kmax + (j-jmin-q-1)*kmax + k-kmin]>0)
        ggcmem[(i-imin)*jmax*kmax + (j-jmin-q-1)*kmax + k-kmin]+=1;

		if(p->gcb4[n][3]==2)
		for(q=0;q<p->margin;++q)
		if(p->flag4[(i-imin)*jmax*kmax + (j-jmin+q+1)*kmax + k-kmin]>0)
		ggcmem[(i-imin)*jmax*kmax + (j-jmin+q+1)*kmax + k-kmin]+=1;

		if(p->gcb4[n][3]==5)
		for(q=0;q<p->margin;++q)
		if(p->flag4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-q-1]>0)
		ggcmem[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-q-1]+=1;

		if(p->gcb4[n][3]==6)
		for(q=0;q<p->margin;++q)
		if(p->flag4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+q+1]>0)
		ggcmem[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+q+1]+=1;
	}

// count entries
	ggccount=0;
	a=0;
	for(i=0;i<imax;++i)
	for(j=0;j<jmax;++j)
	for(k=0;k<kmax;++k)
	{
	if(ggcmem[a]>1)
	{

	++ggccount;
	}

	++a;
	}
	
	p->Iresize(ggc,ggcsize,ggccount*p->margin, 3,3);
	ggcsize=ggccount*p->margin;

//--------------------------
//WALL2
	n=0;
	QQGC4LOOP
	{
        i=p->gcb4[qq][0];
		j=p->gcb4[qq][1];
		k=p->gcb4[qq][2];

		if(p->gcb4[qq][3]==1)
		for(q=0;q<p->margin;++q)
		if(ggcmem[(i-imin-q-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]>1)
        {
            if(ggcmem[(i-imin-q-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]<10)
            {
            ggcmem[(i-imin-q-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]=n+10;
			 ggc[n][0]=i-q-1;
			 ggc[n][1]=j;
			 ggc[n][2]=k;
			 ++n;
            }

        }

        if(p->gcb4[qq][3]==4)
        for(q=0;q<p->margin;++q)
		if(ggcmem[(i-imin+q+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]>1)
        {
            if(ggcmem[(i-imin+q+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]<10)
            {
            ggcmem[(i-imin+q+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]=n+10;
			ggc[n][0]=i+q+1;
			ggc[n][1]=j;
			ggc[n][2]=k;
			++n;
            }

        }

        if(p->gcb4[qq][3]==3)
        for(q=0;q<p->margin;++q)
		if(ggcmem[(i-imin)*jmax*kmax + (j-jmin-q-1)*kmax + k-kmin]>1)
        {
            if(ggcmem[(i-imin)*jmax*kmax + (j-jmin-q-1)*kmax + k-kmin]<10)
            {
            ggcmem[(i-imin)*jmax*kmax + (j-jmin-q-1)*kmax + k-kmin]=n+10;
			ggc[n][0]=i;
			ggc[n][1]=j-q-1;
			ggc[n][2]=k;
			++n;
            }
        }

		if(p->gcb4[qq][3]==2)
		for(q=0;q<p->margin;++q)
		if(ggcmem[(i-imin)*jmax*kmax + (j-jmin+q+1)*kmax + k-kmin]>1)
        {
            if(ggcmem[(i-imin)*jmax*kmax + (j-jmin+q+1)*kmax + k-kmin]<10)
            {
            ggcmem[(i-imin)*jmax*kmax + (j-jmin+q+1)*kmax + k-kmin]=n+10;
			ggc[n][0]=i;
			ggc[n][1]=j+q+1;
			ggc[n][2]=k;

			++n;
            }
        }

        if(p->gcb4[qq][3]==5)
        for(q=0;q<p->margin;++q)
		if(ggcmem[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-q-1]>1)
        {
            if(ggcmem[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-q-1]<10)
            {
            ggcmem[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-q-1]=n+10;
			ggc[n][0]=i;
			ggc[n][1]=j;
			ggc[n][2]=k-q-1;
			++n;
            }
        }

        if(p->gcb4[qq][3]==6)
        for(q=0;q<p->margin;++q)
		if(ggcmem[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+q+1]>1)
        {
            if(ggcmem[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+q+1]<10)
            {
            ggcmem[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+q+1]=n+10;
			ggc[n][0]=i;
			ggc[n][1]=j;
			ggc[n][2]=k+q+1;
			++n;
            }
        }
	}
	ggccount=n;
	
	//cout<<p->mpirank<<" GGC: "<<ggccount<<endl;

	p->del_Iarray(ggcmem,kmax*jmax*imax);
}


