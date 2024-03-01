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

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"
#include"fieldint6.h"

void ghostcell::column_pt4(lexer* p, fdm* a, fieldint &cval4)
{
	n=0;
    LOOP
	{	
	a->C4.p[n] = cval4(i,j,k);
	
	a->C4.n[n] = cval4(i+1,j,k);
	a->C4.s[n] = cval4(i-1,j,k);
	
	a->C4.w[n] = cval4(i,j+1,k);
	a->C4.e[n] = cval4(i,j-1,k);
	
	a->C4.t[n] = cval4(i,j,k+1);
	a->C4.b[n] = cval4(i,j,k-1);
	++n;
	}
	
	GGC4LOOP
    {
    i=p->gcb4[g][0];
    j=p->gcb4[g][1];
    k=p->gcb4[g][2];

        if(p->gcb4[g][3]==1)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C4.n[n] = cval4(i-q,j,k);
			a->C4.s[n] = cval4(i-q-2,j,k);
			}

        ++n;
        }

        if(p->gcb4[g][3]==2)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C4.e[n] = cval4(i,j+q,k);
			a->C4.w[n] = cval4(i,j+q+2,k);
			}
			
        ++n;
        }

        if(p->gcb4[g][3]==3)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C4.e[n] = cval4(i,j-q-2,k); 
			a->C4.w[n] = cval4(i,j-q,k); 
			}

        ++n;
        }

        if(p->gcb4[g][3]==4)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C4.n[n] = cval4(i+q+2,j,k);
			a->C4.s[n] = cval4(i+q,j,k);
			}

		++n;
        }

        if(p->gcb4[g][3]==5)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C4.b[n] = cval4(i,j,k-q-2);
			a->C4.t[n] = cval4(i,j,k-q);
			}

		++n;
        }

        if(p->gcb4[g][3]==6)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C4.b[n] = cval4(i,j,k+q);
			a->C4.t[n] = cval4(i,j,k+q+2);
			}
		
		++n;
        }
    }
	
	for(g=0;g<p->gcpara1_count;++g)
    {
    i=p->gcpara1[g][0];
    j=p->gcpara1[g][1];
    k=p->gcpara1[g][2];
        
        if(p->gcpara1[g][6]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C4.n[n] = cval4(i-q,j,k);
			a->C4.s[n] = cval4(i-q-2,j,k);
			}
 
		++n;
		}
	}
	
	for(g=0;g<p->gcpara2_count;++g)
    {
    i=p->gcpara2[g][0];
    j=p->gcpara2[g][1];
    k=p->gcpara2[g][2];
        
        if(p->gcpara2[g][6]==1)
		for(q=0;q<margin;++q)   
		{
			if(q<margin-1)
			{
			a->C4.e[n] = cval4(i,j+q,k);
			a->C4.w[n] = cval4(i,j+q+2,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara3_count;++g)
    {
    i=p->gcpara3[g][0];
    j=p->gcpara3[g][1];
    k=p->gcpara3[g][2];
        
        if(p->gcpara3[g][6]==1)
		for(q=0;q<margin;++q)   
		{
			if(q<margin-1)
			{
			a->C4.e[n] = cval4(i,j-q-2,k);
			a->C4.w[n] = cval4(i,j-q,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara4_count;++g)
    {
    i=p->gcpara4[g][0];
    j=p->gcpara4[g][1];
    k=p->gcpara4[g][2];
        
        if(p->gcpara4[g][6]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C4.n[n] = cval4(i+q+2,j,k);
			a->C4.s[n] = cval4(i+q,j,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara5_count;++g)
    {
    i=p->gcpara5[g][0];
    j=p->gcpara5[g][1];
    k=p->gcpara5[g][2];
        
        if(p->gcpara5[g][6]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C4.b[n] = cval4(i,j,k-q-2);
			a->C4.t[n] = cval4(i,j,k-q);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara6_count;++g)
    {
    i=p->gcpara6[g][0];
    j=p->gcpara6[g][1];
    k=p->gcpara6[g][2];
        
        if(p->gcpara6[g][6]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C4.b[n] = cval4(i,j,k+q);
			a->C4.t[n] = cval4(i,j,k+q+2);
			}
		
		++n;
		}
	}
}

void ghostcell::column_pt4a(lexer* p, fdm* a, fieldint &cval4a)
{
	n=0;
    ALOOP
	{	
	a->C4a.p[n] = cval4a(i,j,k);
	
	a->C4a.n[n] = cval4a(i+1,j,k);
	a->C4a.s[n] = cval4a(i-1,j,k);
	
	a->C4a.w[n] = cval4a(i,j+1,k);
	a->C4a.e[n] = cval4a(i,j-1,k);
	
	a->C4a.t[n] = cval4a(i,j,k+1);
	a->C4a.b[n] = cval4a(i,j,k-1);
	++n;
	}
	
	GGC4ALOOP
    {
    i=p->gcb4a[g][0];
    j=p->gcb4a[g][1];
    k=p->gcb4a[g][2];

        if(p->gcb4a[g][3]==1)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C4a.n[n] = cval4a(i-q,j,k);
			a->C4a.s[n] = cval4a(i-q-2,j,k);
			}

        ++n;
        }

        if(p->gcb4a[g][3]==2)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C4a.e[n] = cval4a(i,j+q,k);
			a->C4a.w[n] = cval4a(i,j+q+2,k);
			}
			
        ++n;
        }

        if(p->gcb4a[g][3]==3)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C4a.e[n] = cval4a(i,j-q-2,k);
			a->C4a.w[n] = cval4a(i,j-q,k);
			}

        ++n;
        }

        if(p->gcb4a[g][3]==4)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C4a.n[n] = cval4a(i+q+2,j,k);
			a->C4a.s[n] = cval4a(i+q,j,k);
			}

		++n;
        }

        if(p->gcb4a[g][3]==5)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C4a.b[n] = cval4a(i,j,k-q-2);
			a->C4a.t[n] = cval4a(i,j,k-q);
			}

		++n;
        }

        if(p->gcb4a[g][3]==6)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C4a.b[n] = cval4a(i,j,k+q);
			a->C4a.t[n] = cval4a(i,j,k+q+2);
			}
		
		++n;
        }
    }
	
	
	for(g=0;g<p->gcpara1_count;++g)
    {
    i=p->gcpara1[g][0];
    j=p->gcpara1[g][1];
    k=p->gcpara1[g][2];
	
        if(p->gcpara1[g][7]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C4a.n[n] = cval4a(i-q,j,k);
			a->C4a.s[n] = cval4a(i-q-2,j,k);
			}
 
		++n;
		}
	}
	
	for(g=0;g<p->gcpara2_count;++g)
    {
    i=p->gcpara2[g][0];
    j=p->gcpara2[g][1];
    k=p->gcpara2[g][2];
        
        if(p->gcpara2[g][7]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C4a.e[n] = cval4a(i,j+q,k);
			a->C4a.w[n] = cval4a(i,j+q+2,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara3_count;++g)
    {
    i=p->gcpara3[g][0];
    j=p->gcpara3[g][1];
    k=p->gcpara3[g][2];
        
        if(p->gcpara3[g][7]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C4a.e[n] = cval4a(i,j-q-2,k);
			a->C4a.w[n] = cval4a(i,j-q,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara4_count;++g)
    {
    i=p->gcpara4[g][0];
    j=p->gcpara4[g][1];
    k=p->gcpara4[g][2];
        
        if(p->gcpara4[g][7]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C4a.n[n] = cval4a(i+q+2,j,k);
			a->C4a.s[n] = cval4a(i+q,j,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara5_count;++g)
    {
    i=p->gcpara5[g][0];
    j=p->gcpara5[g][1];
    k=p->gcpara5[g][2];
        
        if(p->gcpara5[g][7]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C4a.b[n] = cval4a(i,j,k-q-2);
			a->C4a.t[n] = cval4a(i,j,k-q);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara6_count;++g)
    {
    i=p->gcpara6[g][0];
    j=p->gcpara6[g][1];
    k=p->gcpara6[g][2];
        
        if(p->gcpara6[g][7]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C4a.b[n] = cval4a(i,j,k+q);
			a->C4a.t[n] = cval4a(i,j,k+q+2);
			}
		
		++n;
		}
	}
}

void ghostcell::column_pt6(lexer* p, fdm* a, fieldint &cval6)
{
	n=0;
    BASELOOP
	{	
	a->C6.p[n] = cval6(i,j,k);
	
	a->C6.n[n] = cval6(i+1,j,k);
	a->C6.s[n] = cval6(i-1,j,k);
	
	a->C6.w[n] = cval6(i,j+1,k);
	a->C6.e[n] = cval6(i,j-1,k);
	
	a->C6.t[n] = cval6(i,j,k+1);
	a->C6.b[n] = cval6(i,j,k-1);
	++n;
	}
  
	//int cellnum1=n;
	GGC6LOOP
    {
    i=p->gcb4[g][0];
    j=p->gcb4[g][1];
    k=p->gcb4[g][2];

        if(fabs(p->gcb4[g][3]==1))
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C6.n[n] = cval6(i-q,j,k);
			a->C6.s[n] = cval6(i-q-2,j,k);
			}

        ++n;
        }

        if(fabs(p->gcb4[g][3]==2))
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C6.e[n] = cval6(i,j+q,k);
			a->C6.w[n] = cval6(i,j+q+2,k);
			}
			
        ++n;
        }

        if(fabs(p->gcb4[g][3]==3))
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C6.e[n] = cval6(i,j-q-2,k); 
			a->C6.w[n] = cval6(i,j-q,k); 
			}

        ++n;
        }

        if(fabs(p->gcb4[g][3]==4))
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C6.n[n] = cval6(i+q+2,j,k);
			a->C6.s[n] = cval6(i+q,j,k);
			}

		++n;
        }

        if(fabs(p->gcb4[g][3]==5))
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C6.b[n] = cval6(i,j,k-q-2);
			a->C6.t[n] = cval6(i,j,k-q);
			}

		++n;
        }

        if(fabs(p->gcb4[g][3]==6))
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			a->C6.b[n] = cval6(i,j,k+q);
			a->C6.t[n] = cval6(i,j,k+q+2);
			}
		
		++n;
        }
    }
	


	for(g=0;g<p->gcpara1_count;++g)
    {
    i=p->gcpara1[g][0];
    j=p->gcpara1[g][1];
    k=p->gcpara1[g][2];
        
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C6.n[n] = cval6(i-q,j,k);
			a->C6.s[n] = cval6(i-q-2,j,k);
			}
 
		++n;
		}
	}
	
	for(g=0;g<p->gcpara2_count;++g)
    {
    i=p->gcpara2[g][0];
    j=p->gcpara2[g][1];
    k=p->gcpara2[g][2];
        
		for(q=0;q<margin;++q)   
		{
			if(q<margin-1)
			{
			a->C6.e[n] = cval6(i,j+q,k);
			a->C6.w[n] = cval6(i,j+q+2,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara3_count;++g)
    {
    i=p->gcpara3[g][0];
    j=p->gcpara3[g][1];
    k=p->gcpara3[g][2];
        
		for(q=0;q<margin;++q)   
		{
			if(q<margin-1)
			{
			a->C6.e[n] = cval6(i,j-q-2,k);
			a->C6.w[n] = cval6(i,j-q,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara4_count;++g)
    {
    i=p->gcpara4[g][0];
    j=p->gcpara4[g][1];
    k=p->gcpara4[g][2];
        
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C6.n[n] = cval6(i+q+2,j,k);
			a->C6.s[n] = cval6(i+q,j,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara5_count;++g)
    {
    i=p->gcpara5[g][0];
    j=p->gcpara5[g][1];
    k=p->gcpara5[g][2];
        
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C6.b[n] = cval6(i,j,k-q-2);
			a->C6.t[n] = cval6(i,j,k-q);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara6_count;++g)
    {
    i=p->gcpara6[g][0];
    j=p->gcpara6[g][1];
    k=p->gcpara6[g][2];
        
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			a->C6.b[n] = cval6(i,j,k+q);
			a->C6.t[n] = cval6(i,j,k+q+2);
			}
		
		++n;
		}
	}
}


