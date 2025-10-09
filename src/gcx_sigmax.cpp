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

void ghostcell::sigmax(lexer* p, fdm* a, double *f)
{
    count=0;
    for(n=0;n<p->gcpara1_count;++n)
    {
    i=p->gcpara1[n][0];
    j=p->gcpara1[n][1];
    k=p->gcpara1[n][2];

        for(q=0;q<paramargin;++q)
        {
        send1[count]=f[(i-p->imin+q)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin];
        ++count;
        }
        
    }

    count=0;
    for(n=0;n<p->gcpara2_count;++n)
    {
    i=p->gcpara2[n][0];
    j=p->gcpara2[n][1];
    k=p->gcpara2[n][2];
        
        for(q=0;q<paramargin;++q)
        {
        send2[count]=f[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-q)*p->kmax + k-p->kmin];
        ++count;
        }
    }

    count=0;
    for(n=0;n<p->gcpara3_count;++n)
    {
    i=p->gcpara3[n][0];
    j=p->gcpara3[n][1];
    k=p->gcpara3[n][2];
        
        for(q=0;q<paramargin;++q)
        {
        send3[count]=f[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+q)*p->kmax + k-p->kmin];
        ++count;
        }
    }

    count=0;
    for(n=0;n<p->gcpara4_count;++n)
    {
    i=p->gcpara4[n][0];
    j=p->gcpara4[n][1];
    k=p->gcpara4[n][2];
        
        for(q=0;q<paramargin;++q)
        {
        send4[count]=f[(i-p->imin-q)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin];
        ++count;
        }
    }

    count=0;
    for(n=0;n<p->gcpara5_count;++n)
    {
    i=p->gcpara5[n][0];
    j=p->gcpara5[n][1];
    k=p->gcpara5[n][2];
        
        for(q=0;q<paramargin;++q)
        {
        send5[count]=f[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+q];
        ++count;
        }
    }

    count=0;
    for(n=0;n<p->gcpara6_count;++n)
    {
    i=p->gcpara6[n][0];
    j=p->gcpara6[n][1];
    k=p->gcpara6[n][2];
        
        for(q=0;q<paramargin;++q)
        {
        send6[count]=f[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-q];
        ++count;
        }
    }


//  Communication

    Sendrecv_double(p->gcpara1_count,p->gcpara2_count,p->gcpara3_count,p->gcpara4_count,p->gcpara5_count,p->gcpara6_count);

//  Unpack

    count=0;
    for(n=0;n<p->gcpara1_count;++n)
    {
    i=p->gcpara1[n][0];
    j=p->gcpara1[n][1];
    k=p->gcpara1[n][2];
        
        for(q=0;q<paramargin;++q)
        {
        f[(i-p->imin-q-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]=recv1[count];
        ++count;
        }
    }

    count=0;
    for(n=0;n<p->gcpara2_count;++n)
    {
    i=p->gcpara2[n][0];
    j=p->gcpara2[n][1];
    k=p->gcpara2[n][2];
        
        for(q=0;q<paramargin;++q)
        {
        f[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+q+1)*p->kmax + k-p->kmin]=recv2[count];
        ++count;
        }
    }

    count=0;
    for(n=0;n<p->gcpara3_count;++n)
    {
    i=p->gcpara3[n][0];
    j=p->gcpara3[n][1];
    k=p->gcpara3[n][2];
        
        for(q=0;q<paramargin;++q)
        {
        f[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-q-1)*p->kmax + k-p->kmin]=recv3[count];
        ++count;
        }
    }

    count=0;
    for(n=0;n<p->gcpara4_count;++n)
    {
    i=p->gcpara4[n][0];
    j=p->gcpara4[n][1];
    k=p->gcpara4[n][2];
        
        for(q=0;q<paramargin;++q)
        {
        f[(i-p->imin+q+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]=recv4[count];
        ++count;
        }
    }

    count=0;
    for(n=0;n<p->gcpara5_count;++n)
    {
    i=p->gcpara5[n][0];
    j=p->gcpara5[n][1];
    k=p->gcpara5[n][2];
        
        for(q=0;q<paramargin;++q)
        {
        f[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-q-1]=recv5[count];
        ++count;
        }
    }

    count=0;
    for(n=0;n<p->gcpara6_count;++n)
    {
    i=p->gcpara6[n][0];
    j=p->gcpara6[n][1];
    k=p->gcpara6[n][2];
        
        for(q=0;q<paramargin;++q)
        {
        f[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+q+1]=recv6[count];
        ++count;
        }
    }

}
