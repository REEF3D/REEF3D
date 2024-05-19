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

#include"grid.h"
#include"lexer.h"
#include"fieldint6.h"
#include"cpt.h"

void grid::column_pt1_assign(lexer* p, fieldint &cval1, cpt &C1)
{
	n=0;
    ULOOP
	{	
	C1.p[n] = cval1(i,j,k);

	C1.n[n] = cval1(i+1,j,k);
	C1.s[n] = cval1(i-1,j,k);
	
	C1.w[n] = cval1(i,j+1,k);
	C1.e[n] = cval1(i,j-1,k);
	
	C1.t[n] = cval1(i,j,k+1);
	C1.b[n] = cval1(i,j,k-1);
	++n;
	}

	
	GGC1LOOP
    {
    i=p->gcb1[g][0];
    j=p->gcb1[g][1];
    k=p->gcb1[g][2];

        if(p->gcb1[g][3]==1)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C1.n[n] = cval1(i-q,j,k);
			C1.s[n] = cval1(i-q-2,j,k);
			}

        ++n;
        }

        if(p->gcb1[g][3]==2)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{	
			C1.e[n] = cval1(i,j+q,k);
			C1.w[n] = cval1(i,j+q+2,k);
			}
			
        ++n;
        }

        if(p->gcb1[g][3]==3)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{	
			C1.e[n] = cval1(i,j-q-2,k);
			C1.w[n] = cval1(i,j-q,k);
			}

        ++n;
        }

        if(p->gcb1[g][3]==4)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
            C1.s[n] = cval1(i+q,j,k);
			 C1.n[n] = cval1(i+q+2,j,k);
			
			}
		++n;
        }

        if(p->gcb1[g][3]==5)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C1.b[n] = cval1(i,j,k-q-2);
			C1.t[n] = cval1(i,j,k-q);
			}

		++n;
        }

        if(p->gcb1[g][3]==6)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C1.b[n] = cval1(i,j,k+q);
			C1.t[n] = cval1(i,j,k+q+2);
			}
		
		++n;
        }
    }


	for(g=0;g<p->gcpara1_count;++g)
    {
    i=p->gcpara1[g][0];
    j=p->gcpara1[g][1];
    k=p->gcpara1[g][2];
        
        if(p->gcpara1[g][3]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			C1.n[n] = cval1(i-q,j,k); 
			C1.s[n] = cval1(i-q-2,j,k); 
			}
		++n;
		}
	}
	
	for(g=0;g<p->gcpara2_count;++g)
    {
    i=p->gcpara2[g][0];
    j=p->gcpara2[g][1];
    k=p->gcpara2[g][2];
        
        if(p->gcpara2[g][3]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			C1.e[n] = cval1(i,j+q,k);
			C1.w[n] = cval1(i,j+q+2,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara3_count;++g)
    {
    i=p->gcpara3[g][0];
    j=p->gcpara3[g][1];
    k=p->gcpara3[g][2];
        
        if(p->gcpara3[g][3]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			C1.e[n] = cval1(i,j-q-2,k);
			C1.w[n] = cval1(i,j-q,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara4_count;++g)
    {
    i=p->gcpara4[g][0];
    j=p->gcpara4[g][1];
    k=p->gcpara4[g][2];
        
        if(p->gcpara4[g][3]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			C1.n[n] = cval1(i+q+2,j,k);
			C1.s[n] = cval1(i+q,j,k); 
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara5_count;++g)
    {
    i=p->gcpara5[g][0];
    j=p->gcpara5[g][1];
    k=p->gcpara5[g][2];
        
        if(p->gcpara5[g][3]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			C1.b[n] = cval1(i,j,k-q-2);
			C1.t[n] = cval1(i,j,k-q);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara6_count;++g)
    {
    i=p->gcpara6[g][0];
    j=p->gcpara6[g][1];
    k=p->gcpara6[g][2];
        
        if(p->gcpara6[g][3]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			C1.b[n] = cval1(i,j,k+q);
			C1.t[n] = cval1(i,j,k+q+2);
			}
		
		++n;
		}
	}
    
    //cout<<p->mpirank<<"Cpt1_n: "<<n<<"  p->veclength: "<<p->veclength<<endl;
}

void grid::column_pt2_assign(lexer* p, fieldint &cval2, cpt &C2)
{
	n=0;
    VLOOP
	{	
	C2.p[n] = cval2(i,j,k);

	C2.n[n] = cval2(i+1,j,k);
	C2.s[n] = cval2(i-1,j,k);
	
	C2.w[n] = cval2(i,j+1,k);
	C2.e[n] = cval2(i,j-1,k);
	
	C2.t[n] = cval2(i,j,k+1);
	C2.b[n] = cval2(i,j,k-1);
	++n;
	}
	
	//int cellnum1=n;
	GGC2LOOP
    {
    i=p->gcb2[g][0];
    j=p->gcb2[g][1];
    k=p->gcb2[g][2];

        if(p->gcb2[g][3]==1)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C2.s[n] = cval2(i-q-2,j,k);
			C2.n[n] = cval2(i-q,j,k);
			}

        ++n;
        }

        if(p->gcb2[g][3]==2)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C2.e[n] = cval2(i,j+q,k);
			C2.w[n] = cval2(i,j+q+2,k);
			}
			
        ++n;
        }

        if(p->gcb2[g][3]==3)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C2.e[n] = cval2(i,j-q-2,k);
			C2.w[n] = cval2(i,j-q,k);
			}

        ++n;
        }

        if(p->gcb2[g][3]==4)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C2.s[n] = cval2(i+q,j,k);
			C2.n[n] = cval2(i+q+2,j,k);
			}

		++n;
        }

        if(p->gcb2[g][3]==5)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C2.b[n] = cval2(i,j,k-q-2);
			C2.t[n] = cval2(i,j,k-q);
			}

		++n;
        }

        if(p->gcb2[g][3]==6)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C2.b[n] = cval2(i,j,k+q);
			C2.t[n] = cval2(i,j,k+q+2);
			}
		
		++n;
        }
    }
	

	for(g=0;g<p->gcpara1_count;++g)
    {
    i=p->gcpara1[g][0];
    j=p->gcpara1[g][1];
    k=p->gcpara1[g][2];
        
        if(p->gcpara1[g][4]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			C2.s[n] = cval2(i-q-2,j,k);
			C2.n[n] = cval2(i-q,j,k);
			}
		++n;
		}
	}
	
	for(g=0;g<p->gcpara2_count;++g)
    {
    i=p->gcpara2[g][0];
    j=p->gcpara2[g][1];
    k=p->gcpara2[g][2];
        
        if(p->gcpara2[g][4]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			C2.e[n] = cval2(i,j+q,k);
			C2.w[n] = cval2(i,j+q+2,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara3_count;++g)
    {
    i=p->gcpara3[g][0];
    j=p->gcpara3[g][1];
    k=p->gcpara3[g][2];
        
        if(p->gcpara3[g][4]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			C2.e[n] = cval2(i,j-q-2,k);
			C2.w[n] = cval2(i,j-q,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara4_count;++g)
    {
    i=p->gcpara4[g][0];
    j=p->gcpara4[g][1];
    k=p->gcpara4[g][2];
	
        if(p->gcpara4[g][4]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			C2.s[n] = cval2(i+q,j,k);
			C2.n[n] = cval2(i+q+2,j,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara5_count;++g)
    {
    i=p->gcpara5[g][0];
    j=p->gcpara5[g][1];
    k=p->gcpara5[g][2];
        
        if(p->gcpara5[g][4]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			C2.b[n] = cval2(i,j,k-q-2);
			C2.t[n] = cval2(i,j,k-q);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara6_count;++g)
    {
    i=p->gcpara6[g][0];
    j=p->gcpara6[g][1];
    k=p->gcpara6[g][2];
        
        if(p->gcpara6[g][4]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{	
			C2.b[n] = cval2(i,j,k+q);
			C2.t[n] = cval2(i,j,k+q+2);
			}
		
		++n;
		}
	}
	
}

void grid::column_pt3_assign(lexer* p, fieldint &cval3, cpt &C3)
{
	n=0;
    WLOOP
	{	
	C3.p[n] = cval3(i,j,k);

	C3.n[n] = cval3(i+1,j,k);
	C3.s[n] = cval3(i-1,j,k);
	
	C3.w[n] = cval3(i,j+1,k);
	C3.e[n] = cval3(i,j-1,k);
	
	C3.t[n] = cval3(i,j,k+1);
	C3.b[n] = cval3(i,j,k-1);
	++n;
	}
	
	GGC3LOOP
    {
    i=p->gcb3[g][0];
    j=p->gcb3[g][1];
    k=p->gcb3[g][2];

        if(p->gcb3[g][3]==1)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C3.n[n] = cval3(i-q,j,k);
			C3.s[n] = cval3(i-q-2,j,k);
			}

        ++n;
        }

        if(p->gcb3[g][3]==2)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C3.e[n] = cval3(i,j+q,k);
			C3.w[n] = cval3(i,j+q+2,k);
			}
			
        ++n;
        }

        if(p->gcb3[g][3]==3)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C3.e[n] = cval3(i,j-q-2,k);
			C3.w[n] = cval3(i,j-q,k);
			}

        ++n;
        }

        if(p->gcb3[g][3]==4)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C3.n[n] = cval3(i+q+2,j,k);
			C3.s[n] = cval3(i+q,j,k);
			}

		++n;
        }

        if(p->gcb3[g][3]==5)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C3.t[n] = cval3(i,j,k-q);
			C3.b[n] = cval3(i,j,k-q-2);
			}

		++n;
        }

        if(p->gcb3[g][3]==6)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C3.t[n] = cval3(i,j,k+q+2);
			C3.b[n] = cval3(i,j,k+q);
			}
		
		++n;
        }
    }
	
	for(g=0;g<p->gcpara1_count;++g)
    {
    i=p->gcpara1[g][0];
    j=p->gcpara1[g][1];
    k=p->gcpara1[g][2];
        
        if(p->gcpara1[g][5]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			C3.n[n] = cval3(i-q,j,k);
			C3.s[n] = cval3(i-q-2,j,k);
			}
 
		++n;
		}
	}
	
	for(g=0;g<p->gcpara2_count;++g)
    {
    i=p->gcpara2[g][0];
    j=p->gcpara2[g][1];
    k=p->gcpara2[g][2];
        
        if(p->gcpara2[g][5]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			C3.e[n] = cval3(i,j+q,k);
			C3.w[n] = cval3(i,j+q+2,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara3_count;++g)
    {
    i=p->gcpara3[g][0];
    j=p->gcpara3[g][1];
    k=p->gcpara3[g][2];
        
        if(p->gcpara3[g][5]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			C3.e[n] = cval3(i,j-q-2,k);
			C3.w[n] = cval3(i,j-q,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara4_count;++g)
    {
    i=p->gcpara4[g][0];
    j=p->gcpara4[g][1];
    k=p->gcpara4[g][2];
        
        if(p->gcpara4[g][5]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			C3.n[n] = cval3(i+q+2,j,k);
			C3.s[n] = cval3(i+q,j,k);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara5_count;++g)
    {
    i=p->gcpara5[g][0];
    j=p->gcpara5[g][1];
    k=p->gcpara5[g][2];
        
        if(p->gcpara5[g][5]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			C3.t[n] = cval3(i,j,k-q);
			C3.b[n] = cval3(i,j,k-q-2);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcpara6_count;++g)
    {
    i=p->gcpara6[g][0];
    j=p->gcpara6[g][1];
    k=p->gcpara6[g][2];
        
        if(p->gcpara6[g][5]==1)
		for(q=0;q<margin;++q)   
		{
			if(q<margin-1)
			{
			C3.t[n] = cval3(i,j,k+q+2);
			C3.b[n] = cval3(i,j,k+q);
			}
		
		++n;
		}
	}
}

void grid::column_pt4_assign(lexer* p, fieldint &cval4, cpt &C4)
{
	n=0;
    LOOP
	{	
	C4.p[n] = cval4(i,j,k);
	
	C4.n[n] = cval4(i+1,j,k);
	C4.s[n] = cval4(i-1,j,k);
	
	C4.w[n] = cval4(i,j+1,k);
	C4.e[n] = cval4(i,j-1,k);
	
	C4.t[n] = cval4(i,j,k+1);
	C4.b[n] = cval4(i,j,k-1);
	++n;
	}
	
	//int cellnum1=n;
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
			C4.n[n] = cval4(i-q,j,k);
			C4.s[n] = cval4(i-q-2,j,k);
			}

        ++n;
        }

        if(p->gcb4[g][3]==2)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C4.e[n] = cval4(i,j+q,k);
			C4.w[n] = cval4(i,j+q+2,k);
			}
			
        ++n;
        }

        if(p->gcb4[g][3]==3)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C4.e[n] = cval4(i,j-q-2,k); 
			C4.w[n] = cval4(i,j-q,k); 
			}

        ++n;
        }

        if(p->gcb4[g][3]==4)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C4.n[n] = cval4(i+q+2,j,k);
			C4.s[n] = cval4(i+q,j,k);
			}

		++n;
        }

        if(p->gcb4[g][3]==5)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C4.b[n] = cval4(i,j,k-q-2);
			C4.t[n] = cval4(i,j,k-q);
			}

		++n;
        }

        if(p->gcb4[g][3]==6)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C4.b[n] = cval4(i,j,k+q);
			C4.t[n] = cval4(i,j,k+q+2);
			}
		
		++n;
        }
    }
	
	//cout<<p->mpirank<<" gcbcount4: "<<p->gcb4_count*p->margin<<" N4: "<<n-cellnum1<<endl;
	
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
			C4.n[n] = cval4(i-q,j,k);
			C4.s[n] = cval4(i-q-2,j,k);
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
			C4.e[n] = cval4(i,j+q,k);
			C4.w[n] = cval4(i,j+q+2,k);
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
			C4.e[n] = cval4(i,j-q-2,k);
			C4.w[n] = cval4(i,j-q,k);
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
			C4.n[n] = cval4(i+q+2,j,k);
			C4.s[n] = cval4(i+q,j,k);
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
			C4.b[n] = cval4(i,j,k-q-2);
			C4.t[n] = cval4(i,j,k-q);
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
			C4.b[n] = cval4(i,j,k+q);
			C4.t[n] = cval4(i,j,k+q+2);
			}
		
		++n;
		}
	}
}

void grid::column_pt4a_assign(lexer* p, fieldint &cval4a, cpt &C4a)
{
	n=0;
    ALOOP
	{	
	C4a.p[n] = cval4a(i,j,k);
	
	C4a.n[n] = cval4a(i+1,j,k);
	C4a.s[n] = cval4a(i-1,j,k);
	
	C4a.w[n] = cval4a(i,j+1,k);
	C4a.e[n] = cval4a(i,j-1,k);
	
	C4a.t[n] = cval4a(i,j,k+1);
	C4a.b[n] = cval4a(i,j,k-1);
	++n;
	}
	
	//int cellnum1=n;
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
			C4a.n[n] = cval4a(i-q,j,k);
			C4a.s[n] = cval4a(i-q-2,j,k);
			}

        ++n;
        }

        if(p->gcb4a[g][3]==2)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C4a.e[n] = cval4a(i,j+q,k);
			C4a.w[n] = cval4a(i,j+q+2,k);
			}
			
        ++n;
        }

        if(p->gcb4a[g][3]==3)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C4a.e[n] = cval4a(i,j-q-2,k);
			C4a.w[n] = cval4a(i,j-q,k);
			}

        ++n;
        }

        if(p->gcb4a[g][3]==4)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C4a.n[n] = cval4a(i+q+2,j,k);
			C4a.s[n] = cval4a(i+q,j,k);
			}

		++n;
        }

        if(p->gcb4a[g][3]==5)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C4a.b[n] = cval4a(i,j,k-q-2);
			C4a.t[n] = cval4a(i,j,k-q);
			}

		++n;
        }

        if(p->gcb4a[g][3]==6)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C4a.b[n] = cval4a(i,j,k+q);
			C4a.t[n] = cval4a(i,j,k+q+2);
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
			C4a.n[n] = cval4a(i-q,j,k);
			C4a.s[n] = cval4a(i-q-2,j,k);
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
			C4a.e[n] = cval4a(i,j+q,k);
			C4a.w[n] = cval4a(i,j+q+2,k);
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
			C4a.e[n] = cval4a(i,j-q-2,k);
			C4a.w[n] = cval4a(i,j-q,k);
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
			C4a.n[n] = cval4a(i+q+2,j,k);
			C4a.s[n] = cval4a(i+q,j,k);
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
			C4a.b[n] = cval4a(i,j,k-q-2);
			C4a.t[n] = cval4a(i,j,k-q);
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
			C4a.b[n] = cval4a(i,j,k+q);
			C4a.t[n] = cval4a(i,j,k+q+2);
			}
		
		++n;
		}
	}
}

void grid::column_pt6_assign(lexer* p, fieldint &cval6, cpt &C6)
{
	n=0;
    BASELOOP
	{	
	C6.p[n] = cval6(i,j,k);
	
	C6.n[n] = cval6(i+1,j,k);
	C6.s[n] = cval6(i-1,j,k);
	
	C6.w[n] = cval6(i,j+1,k);
	C6.e[n] = cval6(i,j-1,k);
	
	C6.t[n] = cval6(i,j,k+1);
	C6.b[n] = cval6(i,j,k-1);
	++n;
	}
  
	//int cellnum1=n;
	GGC6LOOP
    {
    i=p->gcb4[g][0];
    j=p->gcb4[g][1];
    k=p->gcb4[g][2];

        if(p->gcb4[g][3]==1)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C6.n[n] = cval6(i-q,j,k);
			C6.s[n] = cval6(i-q-2,j,k);
			}

        ++n;
        }

        if(p->gcb4[g][3]==2)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C6.e[n] = cval6(i,j+q,k);
			C6.w[n] = cval6(i,j+q+2,k);
			}
			
        ++n;
        }

        if(p->gcb4[g][3]==3)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C6.e[n] = cval6(i,j-q-2,k); 
			C6.w[n] = cval6(i,j-q,k); 
			}

        ++n;
        }

        if(p->gcb4[g][3]==4)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C6.n[n] = cval6(i+q+2,j,k);
			C6.s[n] = cval6(i+q,j,k);
			}

		++n;
        }

        if(p->gcb4[g][3]==5)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C6.b[n] = cval6(i,j,k-q-2);
			C6.t[n] = cval6(i,j,k-q);
			}

		++n;
        }

        if(p->gcb4[g][3]==6)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C6.b[n] = cval6(i,j,k+q);
			C6.t[n] = cval6(i,j,k+q+2);
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
			C6.n[n] = cval6(i-q,j,k);
			C6.s[n] = cval6(i-q-2,j,k);
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
			C6.e[n] = cval6(i,j+q,k);
			C6.w[n] = cval6(i,j+q+2,k);
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
			C6.e[n] = cval6(i,j-q-2,k);
			C6.w[n] = cval6(i,j-q,k);
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
			C6.n[n] = cval6(i+q+2,j,k);
			C6.s[n] = cval6(i+q,j,k);
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
			C6.b[n] = cval6(i,j,k-q-2);
			C6.t[n] = cval6(i,j,k-q);
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
			C6.b[n] = cval6(i,j,k+q);
			C6.t[n] = cval6(i,j,k+q+2);
			}
		
		++n;
		}
	}
}

