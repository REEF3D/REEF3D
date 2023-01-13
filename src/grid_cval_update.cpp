/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
#include"fieldint.h"

void grid::cval_update1(lexer* p, fieldint &cval1)
{
    count=0;

    ULOOP
	{
    cval1(i,j,k)=count;
    ++count;
	}

	GC1LOOP
    {
    i=p->gcb1[n][0];
    j=p->gcb1[n][1];
    k=p->gcb1[n][2];
    
        if(p->gcb1[n][3]==1)
        for(q=0;q<margin;++q)
        {
        cval1(i-1-q,j,k)=count;
        ++count;
        }

        if(p->gcb1[n][3]==2)
        for(q=0;q<margin;++q)
        {
        cval1(i,j+1+q,k)=count;
        ++count;
        }

        if(p->gcb1[n][3]==3)
        for(q=0;q<margin;++q)
        {
        cval1(i,j-1-q,k)=count;
        ++count;
        }

        if(p->gcb1[n][3]==4)
        for(q=0;q<margin;++q)
        {
        cval1(i+1+q,j,k)=count;
        ++count;
        }

        if(p->gcb1[n][3]==5)
        for(q=0;q<margin;++q)
        {
        cval1(i,j,k-1-q)=count;
        ++count;
        }

        if(p->gcb1[n][3]==6)
        for(q=0;q<margin;++q)
        {
        cval1(i,j,k+1+q)=count;
        ++count;
        }
    }
	

	for(n=0;n<p->gcpara1_count;++n)
    {
    i=p->gcpara1[n][0];
    j=p->gcpara1[n][1];
    k=p->gcpara1[n][2];
        
        if(p->gcpara1[n][3]==1)
		for(q=0;q<margin;++q)
        {
        cval1(i-1-q,j,k)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara2_count;++n)
    {
    i=p->gcpara2[n][0];
    j=p->gcpara2[n][1];
    k=p->gcpara2[n][2];
        
        if(p->gcpara2[n][3]==1)
		for(q=0;q<margin;++q)
        {
        cval1(i,j+1+q,k)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara3_count;++n)
    {
    i=p->gcpara3[n][0];
    j=p->gcpara3[n][1];
    k=p->gcpara3[n][2];
        
        if(p->gcpara3[n][3]==1)
		for(q=0;q<margin;++q)
        {
        cval1(i,j-1-q,k)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara4_count;++n)
    {
    i=p->gcpara4[n][0];
    j=p->gcpara4[n][1];
    k=p->gcpara4[n][2];
        
        if(p->gcpara4[n][3]==1)
		for(q=0;q<margin;++q)
        {
        cval1(i+1+q,j,k)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara5_count;++n)
    {
    i=p->gcpara5[n][0];
    j=p->gcpara5[n][1];
    k=p->gcpara5[n][2];
        
        if(p->gcpara5[n][3]==1)
		for(q=0;q<margin;++q)
        {
        cval1(i,j,k-1-q)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara6_count;++n)
    {
    i=p->gcpara6[n][0];
    j=p->gcpara6[n][1];
    k=p->gcpara6[n][2];
        
        if(p->gcpara6[n][3]==1)
		for(q=0;q<margin;++q)
        {
        cval1(i,j,k+1+q)=count;
        ++count;
        }
	}
}

void grid::cval_update2(lexer* p, fieldint &cval2)
{
    count=0;

    VLOOP
	{
    cval2(i,j,k)=count;
    ++count;
	}

	GC2LOOP
    {
    i=p->gcb2[n][0];
    j=p->gcb2[n][1];
    k=p->gcb2[n][2];

        if(p->gcb2[n][3]==1)
        for(q=0;q<margin;++q)
        {
        cval2(i-1-q,j,k)=count;
        ++count;
        }

        if(p->gcb2[n][3]==2)
        for(q=0;q<margin;++q)
        {
        cval2(i,j+1+q,k)=count;
        ++count;
        }

        if(p->gcb2[n][3]==3)
        for(q=0;q<margin;++q)
        {
        cval2(i,j-1-q,k)=count;
        ++count;
        }

        if(p->gcb2[n][3]==4)
        for(q=0;q<margin;++q)
        {
        cval2(i+1+q,j,k)=count;
        ++count;
        }

        if(p->gcb2[n][3]==5)
        for(q=0;q<margin;++q)
        {
        cval2(i,j,k-1-q)=count;
        ++count;
        }

        if(p->gcb2[n][3]==6)
        for(q=0;q<margin;++q)
        {
        cval2(i,j,k+1+q)=count;
        ++count;
        }
    }
	
	
	for(n=0;n<p->gcpara1_count;++n)
    {
    i=p->gcpara1[n][0];
    j=p->gcpara1[n][1];
    k=p->gcpara1[n][2];
        
        if(p->gcpara1[n][4]==1)
		for(q=0;q<margin;++q)		
        {
        cval2(i-1-q,j,k)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara2_count;++n)
    {
    i=p->gcpara2[n][0];
    j=p->gcpara2[n][1];
    k=p->gcpara2[n][2];
        
        if(p->gcpara2[n][4]==1)
		for(q=0;q<margin;++q)
        {
        cval2(i,j+1+q,k)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara3_count;++n)
    {
    i=p->gcpara3[n][0];
    j=p->gcpara3[n][1];
    k=p->gcpara3[n][2];
        
        if(p->gcpara3[n][4]==1)
		for(q=0;q<margin;++q)		
        {
        cval2(i,j-1-q,k)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara4_count;++n)
    {
    i=p->gcpara4[n][0];
    j=p->gcpara4[n][1];
    k=p->gcpara4[n][2];
        
        if(p->gcpara4[n][4]==1)
		for(q=0;q<margin;++q)
        {
        cval2(i+1+q,j,k)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara5_count;++n)
    {
    i=p->gcpara5[n][0];
    j=p->gcpara5[n][1];
    k=p->gcpara5[n][2];
        
        if(p->gcpara5[n][4]==1)
		for(q=0;q<margin;++q)		
        {
        cval2(i,j,k-1-q)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara6_count;++n)
    {
    i=p->gcpara6[n][0];
    j=p->gcpara6[n][1];
    k=p->gcpara6[n][2];
        
        if(p->gcpara6[n][4]==1)
		for(q=0;q<margin;++q)
        {
        cval2(i,j,k+1+q)=count;
        ++count;
        }
	}
}

void grid::cval_update3(lexer* p, fieldint &cval3)
{
    count=0;

    WLOOP
	{
    cval3(i,j,k)=count;
    ++count;
	}

	GC3LOOP
    {
    i=p->gcb3[n][0];
    j=p->gcb3[n][1];
    k=p->gcb3[n][2];

        if(p->gcb3[n][3]==1)
        for(q=0;q<margin;++q)
        {
        cval3(i-1-q,j,k)=count;
        ++count;
        }

        if(p->gcb3[n][3]==2)
        for(q=0;q<margin;++q)
        {
        cval3(i,j+1+q,k)=count;
        ++count;
        }

        if(p->gcb3[n][3]==3)
        for(q=0;q<margin;++q)
        {
        cval3(i,j-1-q,k)=count;
        ++count;
        }

        if(p->gcb3[n][3]==4)
        for(q=0;q<margin;++q)
        {
        cval3(i+1+q,j,k)=count;
        ++count;
        }

        if(p->gcb3[n][3]==5)
        for(q=0;q<margin;++q)
        {
        cval3(i,j,k-1-q)=count;
        ++count;
        }

        if(p->gcb3[n][3]==6)
        for(q=0;q<margin;++q)
        {
        cval3(i,j,k+1+q)=count;
        ++count;
        }
	}
	
	
	for(n=0;n<p->gcpara1_count;++n)
    {
    i=p->gcpara1[n][0];
    j=p->gcpara1[n][1];
    k=p->gcpara1[n][2];
        
        if(p->gcpara1[n][5]==1)
		for(q=0;q<margin;++q)
        {
        cval3(i-1-q,j,k)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara2_count;++n)
    {
    i=p->gcpara2[n][0];
    j=p->gcpara2[n][1];
    k=p->gcpara2[n][2];
        
        if(p->gcpara2[n][5]==1)
		for(q=0;q<margin;++q)   
        {
        cval3(i,j+1+q,k)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara3_count;++n)
    {
    i=p->gcpara3[n][0];
    j=p->gcpara3[n][1];
    k=p->gcpara3[n][2];
        
        if(p->gcpara3[n][5]==1)
        for(q=0;q<margin;++q)
        {
        cval3(i,j-1-q,k)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara4_count;++n)
    {
    i=p->gcpara4[n][0];
    j=p->gcpara4[n][1];
    k=p->gcpara4[n][2];
        
        if(p->gcpara4[n][5]==1)
		for(q=0;q<margin;++q)
        {
        cval3(i+1+q,j,k)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara5_count;++n)
    {
    i=p->gcpara5[n][0];
    j=p->gcpara5[n][1];
    k=p->gcpara5[n][2];
        
        if(p->gcpara5[n][5]==1)
		for(q=0;q<margin;++q)
        {
        cval3(i,j,k-1-q)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara6_count;++n)
    {
    i=p->gcpara6[n][0];
    j=p->gcpara6[n][1];
    k=p->gcpara6[n][2];
        
        if(p->gcpara6[n][5]==1)
		for(q=0;q<margin;++q)
        {
        cval3(i,j,k+1+q)=count;
        ++count;
        }
	}
}

void grid::cval_update4(lexer* p, fieldint &cval4)
{
    count=0;

    FLUIDLOOP
	{
    cval4(i,j,k)=count;
    ++count;
	}
	

	GC4LOOP
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];

        if(p->gcb4[n][3]==1)
        for(q=0;q<margin;++q)
        {
        cval4(i-1-q,j,k)=count;
        ++count;
        }

        if(p->gcb4[n][3]==2)
        for(q=0;q<margin;++q)
        {
        cval4(i,j+1+q,k)=count;
        ++count;
        }

        if(p->gcb4[n][3]==3)
        for(q=0;q<margin;++q)
        {
        cval4(i,j-1-q,k)=count;
        ++count;
        }

        if(p->gcb4[n][3]==4)
        for(q=0;q<margin;++q)
        {
        cval4(i+1+q,j,k)=count;
        ++count;
        }

        if(p->gcb4[n][3]==5)
        for(q=0;q<margin;++q)
        {
        cval4(i,j,k-1-q)=count;
        ++count;
        }

        if(p->gcb4[n][3]==6)
        for(q=0;q<margin;++q)
        {
        cval4(i,j,k+1+q)=count;
        ++count;
        }
    }
	
	
	for(n=0;n<p->gcpara1_count;++n)
    {
    i=p->gcpara1[n][0];
    j=p->gcpara1[n][1];
    k=p->gcpara1[n][2];
        
        if(p->gcpara1[n][6]==1)
		for(q=0;q<margin;++q)
        {
        cval4(i-1-q,j,k)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara2_count;++n)
    {
    i=p->gcpara2[n][0];
    j=p->gcpara2[n][1];
    k=p->gcpara2[n][2];
        
        if(p->gcpara2[n][6]==1)
		for(q=0;q<margin;++q)
        {
        cval4(i,j+1+q,k)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara3_count;++n)
    {
    i=p->gcpara3[n][0];
    j=p->gcpara3[n][1];
    k=p->gcpara3[n][2];
        
        if(p->gcpara3[n][6]==1)
		for(q=0;q<margin;++q)
        {
        cval4(i,j-1-q,k)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara4_count;++n)
    {
    i=p->gcpara4[n][0];
    j=p->gcpara4[n][1];
    k=p->gcpara4[n][2];
        
        if(p->gcpara4[n][6]==1)
		for(q=0;q<margin;++q)
        {
        cval4(i+1+q,j,k)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara5_count;++n)
    {
    i=p->gcpara5[n][0];
    j=p->gcpara5[n][1];
    k=p->gcpara5[n][2];
        
        if(p->gcpara5[n][6]==1)
		for(q=0;q<margin;++q)
        {
        cval4(i,j,k-1-q)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara6_count;++n)
    {
    i=p->gcpara6[n][0];
    j=p->gcpara6[n][1];
    k=p->gcpara6[n][2];
        
        if(p->gcpara6[n][6]==1)
		for(q=0;q<margin;++q)
        {
        cval4(i,j,k+1+q)=count;
        ++count;
        }
	}
}

void grid::cval_update4a(lexer* p, fieldint &cval4a)
{
    count=0;

    ALOOP
	{
    cval4a(i,j,k)=count;
    ++count;
	}
	

	GC4ALOOP
    {
    i=p->gcb4a[n][0];
    j=p->gcb4a[n][1];
    k=p->gcb4a[n][2];

        if(p->gcb4a[n][3]==1)
        for(q=0;q<margin;++q)
        {
        cval4a(i-1-q,j,k)=count;
        ++count;
        }

        if(p->gcb4a[n][3]==2)
        for(q=0;q<margin;++q)
        {
        cval4a(i,j+1+q,k)=count;
        ++count;
        }

        if(p->gcb4a[n][3]==3)
        for(q=0;q<margin;++q)
        {
        cval4a(i,j-1-q,k)=count;
        ++count;
        }

        if(p->gcb4a[n][3]==4)
        for(q=0;q<margin;++q)
        {
        cval4a(i+1+q,j,k)=count;
        ++count;
        }

        if(p->gcb4a[n][3]==5)
        for(q=0;q<margin;++q)
        {
        cval4a(i,j,k-1-q)=count;
        ++count;
        }

        if(p->gcb4a[n][3]==6)
        for(q=0;q<margin;++q)
        {
        cval4a(i,j,k+1+q)=count;
        ++count;
        }
    }
	
	
	for(n=0;n<p->gcpara1_count;++n)
    {
    i=p->gcpara1[n][0];
    j=p->gcpara1[n][1];
    k=p->gcpara1[n][2];
        
        if(p->gcpara1[n][7]==1)
		for(q=0;q<margin;++q)
        {
        cval4a(i-1-q,j,k)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara2_count;++n)
    {
    i=p->gcpara2[n][0];
    j=p->gcpara2[n][1];
    k=p->gcpara2[n][2];
        
        if(p->gcpara2[n][7]==1)
		for(q=0;q<margin;++q)
        {
        cval4a(i,j+1+q,k)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara3_count;++n)
    {
    i=p->gcpara3[n][0];
    j=p->gcpara3[n][1];
    k=p->gcpara3[n][2];
        
        if(p->gcpara3[n][7]==1)
		for(q=0;q<margin;++q)
        {
        cval4a(i,j-1-q,k)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara4_count;++n)
    {
    i=p->gcpara4[n][0];
    j=p->gcpara4[n][1];
    k=p->gcpara4[n][2];
        
        if(p->gcpara4[n][7]==1)
		for(q=0;q<margin;++q)
        {
        cval4a(i+1+q,j,k)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara5_count;++n)
    {
    i=p->gcpara5[n][0];
    j=p->gcpara5[n][1];
    k=p->gcpara5[n][2];
        
        if(p->gcpara5[n][7]==1)
		for(q=0;q<margin;++q)
        {
        cval4a(i,j,k-1-q)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara6_count;++n)
    {
    i=p->gcpara6[n][0];
    j=p->gcpara6[n][1];
    k=p->gcpara6[n][2];
        
        if(p->gcpara6[n][7]==1)
		for(q=0;q<margin;++q)
        {
        cval4a(i,j,k+1+q)=count;
        ++count;
        }
	}
}

void grid::cval_update6(lexer* p, fieldint &cval6)
{
    count=0;

    BASELOOP
	{
    cval6(i,j,k)=count;
    ++count;
	}
	

	GC6LOOP
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];

        if(p->gcb4[n][3]==1)
        for(q=0;q<margin;++q)
        {
        cval6(i-1-q,j,k)=count;
        ++count;
        }

        if(p->gcb4[n][3]==2)
        for(q=0;q<margin;++q)
        {
        cval6(i,j+1+q,k)=count;
        ++count;
        }

        if(p->gcb4[n][3]==3)
        for(q=0;q<margin;++q)
        {
        cval6(i,j-1-q,k)=count;
        ++count;
        }

        if(p->gcb4[n][3]==4)
        for(q=0;q<margin;++q)
        {
        cval6(i+1+q,j,k)=count;
        ++count;
        }

        if(p->gcb4[n][3]==5)
        for(q=0;q<margin;++q)
        {
        cval6(i,j,k-1-q)=count;
        ++count;
        }

        if(p->gcb4[n][3]==6)
        for(q=0;q<margin;++q)
        {
        cval6(i,j,k+1+q)=count;
        ++count;
        }
    }
	
	
	for(n=0;n<p->gcpara1_count;++n)
    {
    i=p->gcpara1[n][0];
    j=p->gcpara1[n][1];
    k=p->gcpara1[n][2];
        
        //if(p->gcpara1[n][6]==1)
		for(q=0;q<margin;++q)
        {
        cval6(i-1-q,j,k)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara2_count;++n)
    {
    i=p->gcpara2[n][0];
    j=p->gcpara2[n][1];
    k=p->gcpara2[n][2];
        
       // if(p->gcpara2[n][6]==1)
		for(q=0;q<margin;++q)
        {
        cval6(i,j+1+q,k)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara3_count;++n)
    {
    i=p->gcpara3[n][0];
    j=p->gcpara3[n][1];
    k=p->gcpara3[n][2];
        
       // if(p->gcpara3[n][6]==1)
		for(q=0;q<margin;++q)
        {
        cval6(i,j-1-q,k)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara4_count;++n)
    {
    i=p->gcpara4[n][0];
    j=p->gcpara4[n][1];
    k=p->gcpara4[n][2];
        
        //if(p->gcpara4[n][6]==1)
		for(q=0;q<margin;++q)
        {
        cval6(i+1+q,j,k)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara5_count;++n)
    {
    i=p->gcpara5[n][0];
    j=p->gcpara5[n][1];
    k=p->gcpara5[n][2];
        
        //if(p->gcpara5[n][6]==1)
		for(q=0;q<margin;++q)
        {
        cval6(i,j,k-1-q)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcpara6_count;++n)
    {
    i=p->gcpara6[n][0];
    j=p->gcpara6[n][1];
    k=p->gcpara6[n][2];
        
        //if(p->gcpara6[n][6]==1)
		for(q=0;q<margin;++q)
        {
        cval6(i,j,k+1+q)=count;
        ++count;
        }
	}
}
