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
#include"fdm.h"

void ghostcell::gcparacoxV(lexer* p, double *f, int gcv)
{
    
//  FILL SEND
    for(q=0;q<p->gcparaco1_count;++q)
    {
    i=p->gcparaco1[q][0];
    j=p->gcparaco1[q][1];
    k=p->gcparaco1[q][2];
	send1[q]=f[IJK];
    }

    for(q=0;q<p->gcparaco3_count;++q)
    {
    i=p->gcparaco3[q][0];
    j=p->gcparaco3[q][1];
    k=p->gcparaco3[q][2];
	send3[q]=f[IJK];
    }

	for(q=0;q<p->gcparaco5_count;++q)
	{
    i=p->gcparaco5[q][0];
    j=p->gcparaco5[q][1];
    k=p->gcparaco5[q][2];
	send5[q]=f[IJK];
	}

	for(q=0;q<p->gcparaco4_count;++q)
	{
    i=p->gcparaco4[q][0];
    j=p->gcparaco4[q][1];
    k=p->gcparaco4[q][2];
	send4[q]=f[IJK];
	}

	for(q=0;q<p->gcparaco2_count;++q)
	{
    i=p->gcparaco2[q][0];
    j=p->gcparaco2[q][1];
    k=p->gcparaco2[q][2];
	send2[q]=f[IJK];
	}

	for(q=0;q<p->gcparaco6_count;++q)
	{
    i=p->gcparaco6[q][0];
    j=p->gcparaco6[q][1];
    k=p->gcparaco6[q][2];
	send6[q]=f[IJK];
	}


    Sendrecv6_double(p->gcparaco1_count,p->gcparaco2_count,p->gcparaco3_count,p->gcparaco4_count,p->gcparaco5_count,p->gcparaco6_count);

//  FILL RECEIVE
    for(q=0;q<p->gcparaco1_count;++q)
    {
    i=p->gcparaco1[q][0];
    j=p->gcparaco1[q][1];
    k=p->gcparaco1[q][2];
	f[Im1JK]=recv1[q];
    }

	for(q=0;q<p->gcparaco3_count;++q)
	{
    i=p->gcparaco3[q][0];
    j=p->gcparaco3[q][1];
    k=p->gcparaco3[q][2];
	f[IJm1K]=recv3[q];
	}

    for(q=0;q<p->gcparaco5_count;++q)
    {
    i=p->gcparaco5[q][0];
    j=p->gcparaco5[q][1];
    k=p->gcparaco5[q][2];
	f[IJKm1]=recv5[q];
    }

	for(q=0;q<p->gcparaco4_count;++q)
	{
    i=p->gcparaco4[q][0];
    j=p->gcparaco4[q][1];
    k=p->gcparaco4[q][2];
	f[Ip1JK]=recv4[q];
	}

	for(q=0;q<p->gcparaco2_count;++q)
	{
    i=p->gcparaco2[q][0];
    j=p->gcparaco2[q][1];
    k=p->gcparaco2[q][2];
	f[IJp1K]=recv2[q];
	}

	for(q=0;q<p->gcparaco6_count;++q)
	{
    i=p->gcparaco6[q][0];
    j=p->gcparaco6[q][1];
    k=p->gcparaco6[q][2];
   	f[IJKp1]=recv6[q];
	}
}

