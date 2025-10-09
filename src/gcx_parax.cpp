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

void ghostcell::gcparax(lexer* p,field& f,int gcv)
{
    paramargin=margin;

//  FILL SEND
    count=0;
    for(q=0;q<p->gcpara1_count;++q)
    {
    i=p->gcpara1[q][0]-1;
    j=p->gcpara1[q][1];
    k=p->gcpara1[q][2];
        
        if(p->gcpara1[q][2+gcv]>=1)
        for(n=0;n<paramargin;++n)
        {
        send1[count]=f(i+n+1,j,k);
        ++count;
        }
    }

    count=0;
    for(q=0;q<p->gcpara3_count;++q)
    {
    i=p->gcpara3[q][0];
    j=p->gcpara3[q][1]-1;
    k=p->gcpara3[q][2];
        
        if(p->gcpara3[q][2+gcv]>=1)
        for(n=0;n<paramargin;++n)
        {
        send3[count]=f(i,j+n+1,k);
        ++count;
        }
    }

    count=0;
	for(q=0;q<p->gcpara5_count;++q)
	{
    i=p->gcpara5[q][0];
    j=p->gcpara5[q][1];
    k=p->gcpara5[q][2]-1;
        
        if(p->gcpara5[q][2+gcv]>=1)
        for(n=0;n<paramargin;++n)
        {
        send5[count]=f(i,j,k+n+1);
        ++count;
        }
	}

    count=0;
	for(q=0;q<p->gcpara4_count;++q)
	{
    i=p->gcpara4[q][0]+1;
    j=p->gcpara4[q][1];
    k=p->gcpara4[q][2];
        
        if(p->gcpara4[q][2+gcv]>=1)
        for(n=0;n<paramargin;++n)
        {
        send4[count]=f(i-n-1,j,k);
        ++count;
        }
	}

    count=0;
	for(q=0;q<p->gcpara2_count;++q)
	{
    i=p->gcpara2[q][0];
    j=p->gcpara2[q][1]+1;
    k=p->gcpara2[q][2];
        
        if(p->gcpara2[q][2+gcv]>=1)
        for(n=0;n<paramargin;++n)
        {
        send2[count]=f(i,j-n-1,k);
        ++count;
        }
	}

    count=0;
	for(q=0;q<p->gcpara6_count;++q)
	{
	i=p->gcpara6[q][0];
    j=p->gcpara6[q][1];
    k=p->gcpara6[q][2]+1;
        
        if(p->gcpara6[q][2+gcv]>=1)
        for(n=0;n<paramargin;++n)
        {
        send6[count]=f(i,j,k-n-1);
        ++count;
        }
	}
    
    
    Sendrecv6_double(p->gcpara1_count*paramargin,p->gcpara2_count*paramargin,p->gcpara3_count*paramargin,p->gcpara4_count*paramargin,p->gcpara5_count*paramargin,p->gcpara6_count*paramargin);

//  FILL RECEIVE

    count=0;
    for(q=0;q<p->gcpara1_count;++q)
    {
    i=p->gcpara1[q][0];
    j=p->gcpara1[q][1];
    k=p->gcpara1[q][2];
        
        if(p->gcpara1[q][2+gcv]>=1)
        for(n=0;n<paramargin;++n)
        {
            if(p->gcpara1[q][2+gcv]==1)
            f(i-n-1,j,k)=recv1[count];
            ++count;
        }
    }

    count=0;
	for(q=0;q<p->gcpara3_count;++q)
	{
    i=p->gcpara3[q][0];
    j=p->gcpara3[q][1];
    k=p->gcpara3[q][2];
        
        if(p->gcpara3[q][2+gcv]>=1)
        for(n=0;n<paramargin;++n)
        {
            if(p->gcpara3[q][2+gcv]==1)
            f(i,j-n-1,k)=recv3[count];
            ++count;
        }
	}

	count=0;
    for(q=0;q<p->gcpara5_count;++q)
    {
    i=p->gcpara5[q][0];
    j=p->gcpara5[q][1];
    k=p->gcpara5[q][2];
        
        if(p->gcpara5[q][2+gcv]>=1)
        for(n=0;n<paramargin;++n)
        {
            if(p->gcpara5[q][2+gcv]==1)
            f(i,j,k-n-1)=recv5[count];
            ++count;
        }
    }

    count=0;
	for(q=0;q<p->gcpara4_count;++q)
	{
    i=p->gcpara4[q][0];
    j=p->gcpara4[q][1];
    k=p->gcpara4[q][2];
        
        if(p->gcpara4[q][2+gcv]>=1)
        for(n=0;n<paramargin;++n)
        {
            if(p->gcpara4[q][2+gcv]==1)
            f(i+n+1,j,k)=recv4[count];
            ++count;
        }
	}

    count=0;
	for(q=0;q<p->gcpara2_count;++q)
	{
    i=p->gcpara2[q][0];
    j=p->gcpara2[q][1];
    k=p->gcpara2[q][2];
        
        if(p->gcpara2[q][2+gcv]>=1)
        for(n=0;n<paramargin;++n)
        {
            if(p->gcpara2[q][2+gcv]==1)
            f(i,j+n+1,k)=recv2[count];
            ++count;
        }
	}

    count=0;
	for(q=0;q<p->gcpara6_count;++q)
	{
    i=p->gcpara6[q][0];
    j=p->gcpara6[q][1];
    k=p->gcpara6[q][2];
        
        if(p->gcpara6[q][2+gcv]>=1)
        for(n=0;n<paramargin;++n)
        {  
            if(p->gcpara6[q][2+gcv]==1)
            f(i,j,k+n+1)=recv6[count];
            ++count;
        }
	}

}

