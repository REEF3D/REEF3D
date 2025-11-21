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

void ghostcell::gcslparax_fh(lexer* p, slice& f, int gcv)
{
    paramargin=margin;

    //  FILL SEND
    count=0;
    for(q=0;q<p->gcslpara1_count;++q)
    {
        i=p->gcslpara1[q][0];
        j=p->gcslpara1[q][1];

        for(n=0;n<paramargin;++n)
        {
            send1[count]=f(i-n-1,j);
            ++count;
        }
    }

    count=0;
    for(q=0;q<p->gcslpara3_count;++q)
    {
        i=p->gcslpara3[q][0];
        j=p->gcslpara3[q][1];

        for(n=0;n<paramargin;++n)
        {
            send3[count]=f(i,j-n-1);
            ++count;
        }
    }

    count=0;
    for(q=0;q<p->gcslpara4_count;++q)
    {
        i=p->gcslpara4[q][0];
        j=p->gcslpara4[q][1];

        for(n=0;n<paramargin;++n)
        {
            send4[count]=f(i+n+1,j);
            ++count;
        }
    }

    count=0;
    for(q=0;q<p->gcslpara2_count;++q)
    {
        i=p->gcslpara2[q][0];
        j=p->gcslpara2[q][1];

        for(n=0;n<paramargin;++n)
        {
            send2[count]=f(i,j+n+1);
            ++count;
        }
    }

    Sendrecv_double(p->gcslpara1_count*paramargin,p->gcslpara2_count*paramargin,p->gcslpara3_count*paramargin,p->gcslpara4_count*paramargin,0,0);

    //  FILL RECEIVE
    count=0;
    for(q=0;q<p->gcslpara1_count;++q)
    {
        i=p->gcslpara1[q][0];
        j=p->gcslpara1[q][1];

        for(n=0;n<paramargin;++n)
        {
            f(i+n,j)+=recv1[count];
            ++count;
        }
    }

    count=0;
    for(q=0;q<p->gcslpara3_count;++q)
    {
        i=p->gcslpara3[q][0];
        j=p->gcslpara3[q][1];

        for(n=0;n<paramargin;++n)
        {
            f(i,j+n)+=recv3[count];
            ++count;
        }
    }

    count=0;
    for(q=0;q<p->gcslpara4_count;++q)
    {
        i=p->gcslpara4[q][0];
        j=p->gcslpara4[q][1];

        for(n=0;n<paramargin;++n)
        {
            f(i-n,j)+=recv4[count];
            ++count;
        }
    }

    count=0;
    for(q=0;q<p->gcslpara2_count;++q)
    {
        i=p->gcslpara2[q][0];
        j=p->gcslpara2[q][1];

        for(n=0;n<paramargin;++n)
        {
            f(i,j-n)+=recv2[count];
            ++count;
        }
    }
}
