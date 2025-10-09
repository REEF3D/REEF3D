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

void ghostcell::gcparaxijk(lexer* p, double *f, int gcv)
{
	starttime=timer();
	
    paramargin=3;
	
//  FILL SEND
    count=0;
    for(q=0;q<p->gcpara1_count;++q)
    {
    i=p->gcpara1[q][0];
    j=p->gcpara1[q][1];
    k=p->gcpara1[q][2];
    
        send1[count]=f[IJK];
        ++count;

        send1[count]=f[Ip1JK];
        ++count;
        
        send1[count]=f[Ip2JK];
        ++count;
    }
	
	count=0;
	for(q=0;q<p->gcpara2_count;++q)
	{
    i=p->gcpara2[q][0];
    j=p->gcpara2[q][1];
    k=p->gcpara2[q][2];
    
        send2[count]=f[IJK];
        ++count;

        send2[count]=f[IJm1K];
        ++count;
  
        send2[count]=f[IJm2K];
        ++count;
	}

    count=0;
    for(q=0;q<p->gcpara3_count;++q)
    {
    i=p->gcpara3[q][0];
    j=p->gcpara3[q][1];
    k=p->gcpara3[q][2];
    
        send3[count]=f[IJK];
        ++count;
        
        send3[count]=f[IJp1K];
        ++count;
     
        send3[count]=f[IJp2K];
        ++count;
    }
	
	count=0;
	for(q=0;q<p->gcpara4_count;++q)
	{
    i=p->gcpara4[q][0];
    j=p->gcpara4[q][1];
    k=p->gcpara4[q][2];
    
        send4[count]=f[IJK];
        ++count;

        send4[count]=f[Im1JK];
        ++count;

        send4[count]=f[Im2JK];
        ++count;
	}

    count=0;
	for(q=0;q<p->gcpara5_count;++q)
	{
    i=p->gcpara5[q][0];
    j=p->gcpara5[q][1];
    k=p->gcpara5[q][2];
    
        send5[count]=f[IJK];
        ++count;

        send5[count]=f[IJKp1];
        ++count;

        send5[count]=f[IJKp2];
        ++count;
	}

    count=0;
	for(q=0;q<p->gcpara6_count;++q)
	{
    i=p->gcpara6[q][0];
    j=p->gcpara6[q][1];
    k=p->gcpara6[q][2];
    
        send6[count]=f[IJK];
        ++count;

        send6[count]=f[IJKm1];
        ++count;

        send6[count]=f[IJKm2];
        ++count;
	}


    Sendrecv6_double(p->gcpara1_count*paramargin,p->gcpara2_count*paramargin,p->gcpara3_count*paramargin,p->gcpara4_count*paramargin,p->gcpara5_count*paramargin,p->gcpara6_count*paramargin);

//  FILL RECEIVE

    count=0;
    for(q=0;q<p->gcpara1_count;++q)
    {
        f[Im1JK]=recv1[count];
        ++count;

        f[Im2JK]=recv1[count];
        ++count;

        f[Im3JK]=recv1[count];
        ++count; 
    }

    count=0;
	for(q=0;q<p->gcpara2_count;++q)
	{
        f[IJp1K]=recv2[count];
        ++count;

        f[IJp2K]=recv2[count];
        ++count;
        
        f[IJp3K]=recv2[count];
        ++count;
	}	
	
	count=0;
	for(q=0;q<p->gcpara3_count;++q)
	{
        f[IJm1K]=recv3[count];
        ++count;

        f[IJm2K]=recv3[count];
        ++count;
        
        f[IJm3K]=recv3[count];
        ++count;
	}

    count=0;
	for(q=0;q<p->gcpara4_count;++q)
	{
        f[Ip1JK]=recv4[count];
        ++count;

        f[Ip2JK]=recv4[count];
        ++count;
        
        f[Ip3JK]=recv4[count];
        ++count;
	}
	
	count=0;
    for(q=0;q<p->gcpara5_count;++q)
    {
        f[IJKm1]=recv5[count];
        ++count;
    
        f[IJKm2]=recv5[count];
        ++count;

        f[IJKm3]=recv5[count];
        ++count;
    }

    count=0;
	for(q=0;q<p->gcpara6_count;++q)
	{
        f[IJKp1]=recv6[count];
        ++count;
  
        f[IJKp2]=recv6[count];
        ++count;

        f[IJKp3]=recv6[count];
        ++count;
	}
	endtime=timer();
	p->xtime+=endtime-starttime;
}

