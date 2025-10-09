/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"ghostcell.h"
#include"lexer.h"
#include"sliceint.h"

void ghostcell::gcslparacox_int(lexer* p,sliceint& f,int gcv)
{
//  FILL SEND
    for(q=0;q<p->gcslparaco1_count;++q)
    {
    i=p->gcslparaco1[q][0];
    j=p->gcslparaco1[q][1];
	isend1[q]=f(i,j);
    }

    for(q=0;q<p->gcslparaco3_count;++q)
    {
    i=p->gcslparaco3[q][0];
    j=p->gcslparaco3[q][1];
	isend3[q]=f(i,j);
    }

	for(q=0;q<p->gcslparaco4_count;++q)
	{
    i=p->gcslparaco4[q][0];
    j=p->gcslparaco4[q][1];
	isend4[q]=f(i,j);
	}

	for(q=0;q<p->gcslparaco2_count;++q)
	{
    i=p->gcslparaco2[q][0];
    j=p->gcslparaco2[q][1];
	isend2[q]=f(i,j);
	}

    Sendrecv_int(p->gcslparaco1_count, p->gcslparaco2_count, p->gcslparaco3_count, p->gcslparaco4_count, 0, 0);

//  FILL RECEIVE
    for(q=0;q<p->gcslparaco1_count;++q)
    {
    i=p->gcslparaco1[q][0];
    j=p->gcslparaco1[q][1];
	f(i-1,j)=irecv1[q];
    }

	for(q=0;q<p->gcslparaco3_count;++q)
	{
    i=p->gcslparaco3[q][0];
    j=p->gcslparaco3[q][1];
	f(i,j-1)=irecv3[q];
	}

	for(q=0;q<p->gcslparaco4_count;++q)
	{
    i=p->gcslparaco4[q][0];
    j=p->gcslparaco4[q][1];
	f(i+1,j)=irecv4[q];
	}

	for(q=0;q<p->gcslparaco2_count;++q)
	{
    i=p->gcslparaco2[q][0];
    j=p->gcslparaco2[q][1];
	f(i,j+1)=irecv2[q];
	}
}

void ghostcell::gcslparacoxV_int(lexer* p, int *f, int gcv)
{
//  FILL SEND
    for(q=0;q<p->gcslparaco1_count;++q)
    {
    i=p->gcslparaco1[q][0];
    j=p->gcslparaco1[q][1];
	isend1[q]=f[IJ];
    }

    for(q=0;q<p->gcslparaco3_count;++q)
    {
    i=p->gcslparaco3[q][0];
    j=p->gcslparaco3[q][1];
	isend3[q]=f[IJ];
    }

	for(q=0;q<p->gcslparaco4_count;++q)
	{
    i=p->gcslparaco4[q][0];
    j=p->gcslparaco4[q][1];
	isend4[q]=f[IJ];
	}

	for(q=0;q<p->gcslparaco2_count;++q)
	{
    i=p->gcslparaco2[q][0];
    j=p->gcslparaco2[q][1];
	isend2[q]=f[IJ];
	}

    Sendrecv_int(p->gcslparaco1_count, p->gcslparaco2_count, p->gcslparaco3_count, p->gcslparaco4_count, 0, 0);

//  FILL RECEIVE
    for(q=0;q<p->gcslparaco1_count;++q)
    {
    i=p->gcslparaco1[q][0];
    j=p->gcslparaco1[q][1];
	f[Im1J]=irecv1[q];
    }

	for(q=0;q<p->gcslparaco3_count;++q)
	{
    i=p->gcslparaco3[q][0];
    j=p->gcslparaco3[q][1];
	f[IJm1]=irecv3[q];
	}

	for(q=0;q<p->gcslparaco4_count;++q)
	{
    i=p->gcslparaco4[q][0];
    j=p->gcslparaco4[q][1];
	f[Ip1J]=irecv4[q];
	}

	for(q=0;q<p->gcslparaco2_count;++q)
	{
    i=p->gcslparaco2[q][0];
    j=p->gcslparaco2[q][1];
	f[IJp1]=irecv2[q];
	}
}

