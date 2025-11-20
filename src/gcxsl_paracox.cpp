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
#include"slice.h"

void ghostcell::gcslparacox(lexer* p, slice& f, int gcv)
{
    //  FILL SEND
    for(q=0;q<p->gcslparaco1_count;++q)
    {
        i=p->gcslparaco1[q][0];
        j=p->gcslparaco1[q][1];
        send1[q]=f(i,j);
    }

    for(q=0;q<p->gcslparaco3_count;++q)
    {
        i=p->gcslparaco3[q][0];
        j=p->gcslparaco3[q][1];
        send3[q]=f(i,j);
    }

    for(q=0;q<p->gcslparaco4_count;++q)
    {
        i=p->gcslparaco4[q][0];
        j=p->gcslparaco4[q][1];
        send4[q]=f(i,j);
    }

    for(q=0;q<p->gcslparaco2_count;++q)
    {
        i=p->gcslparaco2[q][0];
        j=p->gcslparaco2[q][1];
        send2[q]=f(i,j);
    }

    Sendrecv_double(p->gcslparaco1_count,p->gcslparaco2_count,p->gcslparaco3_count,p->gcslparaco4_count,0,0);

    //  FILL RECEIVE
    for(q=0;q<p->gcslparaco1_count;++q)
    {
        i=p->gcslparaco1[q][0];
        j=p->gcslparaco1[q][1];
        f(i-1,j)=recv1[q];
    }

    for(q=0;q<p->gcslparaco3_count;++q)
    {
        i=p->gcslparaco3[q][0];
        j=p->gcslparaco3[q][1];
        f(i,j-1)=recv3[q];
    }

    for(q=0;q<p->gcslparaco4_count;++q)
    {
        i=p->gcslparaco4[q][0];
        j=p->gcslparaco4[q][1];
        f(i+1,j)=recv4[q];
    }

    for(q=0;q<p->gcslparaco2_count;++q)
    {
        i=p->gcslparaco2[q][0];
        j=p->gcslparaco2[q][1];
        f(i,j+1)=recv2[q];
    }
}
