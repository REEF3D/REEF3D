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

#include"mgc2.h"
#include"cart2.h"
#include"lexer.h"

void mgc2::make_ggc(lexer* p)
{
	cart2::ggcsize=1;
	p->Iarray(cart2::ggc,cart2::ggcsize,3);
}

void mgc2::fill_ggc(lexer* p)
{
	int q,qq,n,nn,a;
	int check;
	
	p->Iarray(cart2::ggcmem,kmax*jmax*imax);

//--------------------------
//WALL1

	GC2LOOP
	{
	    i=p->gcb2[n][0];
		j=p->gcb2[n][1];
		k=p->gcb2[n][2];

		if(p->gcb2[n][3]==1)
		for(q=0;q<p->margin;++q)
        if(p->flag2[(i+q+1-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]<0)
        cart2::ggcmem[(i-imin-q-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]+=1;

		if(p->gcb2[n][3]==4)
		for(q=0;q<p->margin;++q)
        if(p->flag2[(i-q-1-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]<0)
		cart2::ggcmem[(i-imin+q+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]+=1;

		if(p->gcb2[n][3]==3)
		for(q=0;q<p->margin;++q)
        if(p->flag2[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-q-1)*p->kmax + k-p->kmin]<0)
        cart2::ggcmem[(i-imin)*jmax*kmax + (j-jmin-q-1)*kmax + k-kmin]+=1;

		if(p->gcb2[n][3]==2)
		for(q=0;q<p->margin;++q)
        if(p->flag2[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+q+1)*p->kmax + k-p->kmin]<0)
		cart2::ggcmem[(i-imin)*jmax*kmax + (j-jmin+q+1)*kmax + k-kmin]+=1;

		if(p->gcb2[n][3]==5)
		for(q=0;q<p->margin;++q)
        if(p->flag2[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-q-1]<0)
		cart2::ggcmem[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-q-1]+=1;

		if(p->gcb2[n][3]==6)
		for(q=0;q<p->margin;++q)
        if(p->flag2[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+q+1]<0)
		cart2::ggcmem[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+q+1]+=1;
	}

// count entries
	cart2::ggccount=0;
	a=0;
	for(i=0;i<imax;++i)
	for(j=0;j<jmax;++j)
	for(k=0;k<kmax;++k)
	{
        if(cart2::ggcmem[a]>1)
        {
        ++cart2::ggccount;
        }

	++a;
	}
	
	p->Iresize(cart2::ggc,cart2::ggcsize,cart2::ggccount*p->margin, 3, 3);
	cart2::ggcsize=cart2::ggccount*p->margin;

//--------------------------
//WALL2
	n=0;
	QQGC2LOOP
	{
        i=p->gcb2[qq][0];
		j=p->gcb2[qq][1];
		k=p->gcb2[qq][2];

		if(p->gcb2[qq][3]==1)
		for(q=0;q<p->margin;++q)
		if(cart2::ggcmem[(i-imin-q-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]>1)
        {
            if(cart2::ggcmem[(i-imin-q-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]<10)
            {
            cart2::ggcmem[(i-imin-q-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]=n+10;
			cart2::ggc[n][0]=i-q-1;
			cart2::ggc[n][1]=j;
			cart2::ggc[n][2]=k;
			++n;
            }
        }

        if(p->gcb2[qq][3]==4)
        for(q=0;q<p->margin;++q)
		if(cart2::ggcmem[(i-imin+q+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]>1)
        {
            if(cart2::ggcmem[(i-imin+q+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]<10)
            {
            cart2::ggcmem[(i-imin+q+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]=n+10;
			cart2::ggc[n][0]=i+q+1;
			cart2::ggc[n][1]=j;
			cart2::ggc[n][2]=k;
			++n;
            }
        }

        if(p->gcb2[qq][3]==3)
        for(q=0;q<p->margin;++q)
		if(cart2::ggcmem[(i-imin)*jmax*kmax + (j-jmin-q-1)*kmax + k-kmin]>1)
        {
            if(cart2::ggcmem[(i-imin)*jmax*kmax + (j-jmin-q-1)*kmax + k-kmin]<10)
            {
            cart2::ggcmem[(i-imin)*jmax*kmax + (j-jmin-q-1)*kmax + k-kmin]=n+10;
			cart2::ggc[n][0]=i;
			cart2::ggc[n][1]=j-q-1;
			cart2::ggc[n][2]=k;
			++n;
            }
        }

		if(p->gcb2[qq][3]==2)
		for(q=0;q<p->margin;++q)
		if(cart2::ggcmem[(i-imin)*jmax*kmax + (j-jmin+q+1)*kmax + k-kmin]>1)
        {
            if(cart2::ggcmem[(i-imin)*jmax*kmax + (j-jmin+q+1)*kmax + k-kmin]<10)
            {
            cart2::ggcmem[(i-imin)*jmax*kmax + (j-jmin+q+1)*kmax + k-kmin]=n+10;
			cart2::ggc[n][0]=i;
			cart2::ggc[n][1]=j+q+1;
			cart2::ggc[n][2]=k;
			++n;
            }

        }

        if(p->gcb2[qq][3]==5)
        for(q=0;q<p->margin;++q)
		if(cart2::ggcmem[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-q-1]>1)
        {
            if(cart2::ggcmem[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-q-1]<10)
            {
            cart2::ggcmem[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-q-1]=n+10;
			cart2::ggc[n][0]=i;
			cart2::ggc[n][1]=j;
			cart2::ggc[n][2]=k-q-1;
			++n;
            }
        }

        if(p->gcb2[qq][3]==6)
        for(q=0;q<p->margin;++q)
		if(cart2::ggcmem[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+q+1]>1)
        {
            if(cart2::ggcmem[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+q+1]<10)
            {
            cart2::ggcmem[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+q+1]=n+10;
			cart2::ggc[n][0]=i;
			cart2::ggc[n][1]=j;
			cart2::ggc[n][2]=k+q+1;
			++n;
            }
        }
	}
	cart2::ggccount=n;


	p->del_Iarray(cart2::ggcmem,kmax*jmax*imax);
}


