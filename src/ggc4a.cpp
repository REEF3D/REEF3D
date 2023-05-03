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

#include"mgc4a.h"
#include"cart4a.h"
#include"lexer.h"

void mgc4a::make_ggc(lexer* p)
{
    cart4a::ggcsize=1;
	p->Iarray(cart4a::ggc,cart4a::ggcsize,3);
}

void mgc4a::fill_ggc(lexer* p)
{
	int q,qq,n,nn,a,check;

	p->Iarray(cart4a::ggcmem,kmax*jmax*imax);
	
//--------------------------
//WALL1

	GC4ALOOP
	{
	    i=p->gcb4a[n][0];
		j=p->gcb4a[n][1];
		k=p->gcb4a[n][2];

		if(fabs(p->gcb4a[n][3])==1)
		for(q=0;q<p->margin;++q)
        cart4a::ggcmem[(i-imin-q-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]+=1;

		if(fabs(p->gcb4a[n][3])==4)
		for(q=0;q<p->margin;++q)
		cart4a::ggcmem[(i-imin+q+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]+=1;

		if(fabs(p->gcb4a[n][3])==3)
		for(q=0;q<p->margin;++q)
        cart4a::ggcmem[(i-imin)*jmax*kmax + (j-jmin-q-1)*kmax + k-kmin]+=1;

		if(fabs(p->gcb4a[n][3])==2)
		for(q=0;q<p->margin;++q)
		cart4a::ggcmem[(i-imin)*jmax*kmax + (j-jmin+q+1)*kmax + k-kmin]+=1;

		if(fabs(p->gcb4a[n][3])==5)
		for(q=0;q<p->margin;++q)
		cart4a::ggcmem[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-q-1]+=1;

		if(fabs(p->gcb4a[n][3])==6)
		for(q=0;q<p->margin;++q)
		cart4a::ggcmem[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+q+1]+=1;
	}

// count entries
	cart4a::ggccount=0;
	a=0;
	for(i=0;i<imax;++i)
	for(j=0;j<jmax;++j)
	for(k=0;k<kmax;++k)
	{
	if(cart4a::ggcmem[a]>1)
	++cart4a::ggccount;

	++a;
	}

	p->Iresize(cart4a::ggc,cart4a::ggcsize,cart4a::ggccount*p->margin, 3, 3);
	cart4a::ggcsize=cart4a::ggccount*p->margin;
	
	

//--------------------------
//WALL2
	n=0;
	QQGC4ALOOP
	{
        i=p->gcb4a[qq][0];
		j=p->gcb4a[qq][1];
		k=p->gcb4a[qq][2];

		if(fabs(p->gcb4a[qq][3])==1)
		for(q=0;q<p->margin;++q)
		if(cart4a::ggcmem[(i-imin-q-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]>1)
        {
            if(cart4a::ggcmem[(i-imin-q-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]<10)
            {
            cart4a::ggcmem[(i-imin-q-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]=n+10;
			cart4a::ggc[n][0]=i-q-1;
			cart4a::ggc[n][1]=j;
			cart4a::ggc[n][2]=k;
			++n;
            }
        }

        if(fabs(p->gcb4a[qq][3])==4)
        for(q=0;q<p->margin;++q)
		if(cart4a::ggcmem[(i-imin+q+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]>1)
        {
            if(cart4a::ggcmem[(i-imin+q+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]<10)
            {
            cart4a::ggcmem[(i-imin+q+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]=n+10;
			cart4a::ggc[n][0]=i+q+1;
			cart4a::ggc[n][1]=j;
			cart4a::ggc[n][2]=k;
			++n;
            }
        }

        if(fabs(p->gcb4a[qq][3])==3)
        for(q=0;q<p->margin;++q)
		if(cart4a::ggcmem[(i-imin)*jmax*kmax + (j-jmin-q-1)*kmax + k-kmin]>1)
        {
            if(cart4a::ggcmem[(i-imin)*jmax*kmax + (j-jmin-q-1)*kmax + k-kmin]<10)
            {
            cart4a::ggcmem[(i-imin)*jmax*kmax + (j-jmin-q-1)*kmax + k-kmin]=n+10;
			cart4a::ggc[n][0]=i;
			cart4a::ggc[n][1]=j-q-1;
			cart4a::ggc[n][2]=k;
			++n;
            }
        }


		if(fabs(p->gcb4a[qq][3])==2)
		for(q=0;q<p->margin;++q)
		if(cart4a::ggcmem[(i-imin)*jmax*kmax + (j-jmin+q+1)*kmax + k-kmin]>1)
        {
            if(cart4a::ggcmem[(i-imin)*jmax*kmax + (j-jmin+q+1)*kmax + k-kmin]<10)
            {
            cart4a::ggcmem[(i-imin)*jmax*kmax + (j-jmin+q+1)*kmax + k-kmin]=n+10;
			cart4a::ggc[n][0]=i;
			cart4a::ggc[n][1]=j+q+1;
			cart4a::ggc[n][2]=k;
			++n;
            }
        }

        if(fabs(p->gcb4a[qq][3])==5)
        for(q=0;q<p->margin;++q)
		if(cart4a::ggcmem[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-q-1]>1)
        {
            if(cart4a::ggcmem[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-q-1]<10)
            {
            cart4a::ggcmem[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-q-1]=n+10;
			cart4a::ggc[n][0]=i;
			cart4a::ggc[n][1]=j;
			cart4a::ggc[n][2]=k-q-1;
			++n;
            }
        }

        if(fabs(p->gcb4a[qq][3])==6)
        for(q=0;q<p->margin;++q)
		if(cart4a::ggcmem[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+q+1]>1)
        {
            if(cart4a::ggcmem[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+q+1]<10)
            {
            cart4a::ggcmem[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+q+1]=n+10;
			cart4a::ggc[n][0]=i;
			cart4a::ggc[n][1]=j;
			cart4a::ggc[n][2]=k+q+1;
			++n;
            }
        }
	}
	cart4a::ggccount=n;
	
	p->del_Iarray(cart4a::ggcmem,kmax*jmax*imax);

}




