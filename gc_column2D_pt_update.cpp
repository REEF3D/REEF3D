/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include"lexer.h"
#include"sliceint1.h"
#include"sliceint2.h"
#include"sliceint4.h"
#include"cpt2D.h"

void ghostcell::column2D_pt1_update(lexer* p, cpt2D &C1)
{
    sliceint1 cval1(p);
    cval2Dupdate1(p,cval1);	
	column2D_pt1(p,cval1,C1);
    cval_gcslpara1(p,cval1);
}

void ghostcell::column2D_pt2_update(lexer* p, cpt2D &C2)
{
    sliceint2 cval2(p);
    cval2Dupdate2(p,cval2);	
	column2D_pt2(p,cval2,C2);
    cval_gcslpara2(p,cval2);
}

void ghostcell::column2D_pt4_update(lexer* p, cpt2D &C4)
{
    sliceint4 cval4(p);
    cval2Dupdate4(p,cval4);	
	column2D_pt4(p,cval4,C4);
    cval_gcslpara4(p,cval4);
}

void ghostcell::column2D_pt1(lexer* p, sliceint &cval1, cpt2D &C1)
{
	n=0;
    SLICELOOP1
	{	
	C1.p[n] = cval1(i,j);
	
	C1.n[n] = cval1(i+1,j);
	C1.s[n] = cval1(i-1,j);
	
	C1.w[n] = cval1(i,j+1);
	C1.e[n] = cval1(i,j-1);

	++n;
	}
	
	GGCSL1LOOP
    {
    i=p->gcbsl1[g][0];
    j=p->gcbsl1[g][1];

        if(p->gcbsl1[g][3]==1)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C1.n[n] = cval1(i-q,j);
			C1.s[n] = cval1(i-q-2,j);
			}

        ++n;
        }

        if(p->gcbsl1[g][3]==2)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C1.e[n] = cval1(i,j+q);
			C1.w[n] = cval1(i,j+q+2);
			}
			
        ++n;
        }

        if(p->gcbsl1[g][3]==3)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C1.e[n] = cval1(i,j-q-2); 
			C1.w[n] = cval1(i,j-q); 
			}

        ++n;
        }

        if(p->gcbsl1[g][3]==4)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C1.n[n] = cval1(i+q+2,j);
			C1.s[n] = cval1(i+q,j);
			}

		++n;
        }
    }
	
	for(g=0;g<p->gcslpara1_count;++g)
    {
    i=p->gcslpara1[g][0];
    j=p->gcslpara1[g][1];
        
        if(p->gcslpara1[g][3]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			C1.n[n] = cval1(i-q,j);
			C1.s[n] = cval1(i-q-2,j);
			}
 
		++n;
		}
	}
	
	for(g=0;g<p->gcslpara2_count;++g)
    {
    i=p->gcslpara2[g][0];
    j=p->gcslpara2[g][1];
        
        if(p->gcslpara2[g][3]==1)
		for(q=0;q<margin;++q)   
		{
			if(q<margin-1)
			{
			C1.e[n] = cval1(i,j+q);
			C1.w[n] = cval1(i,j+q+2);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcslpara3_count;++g)
    {
    i=p->gcslpara3[g][0];
    j=p->gcslpara3[g][1];
        
        if(p->gcslpara3[g][3]==1)
		for(q=0;q<margin;++q)   
		{
			if(q<margin-1)
			{
			C1.e[n] = cval1(i,j-q-2);
			C1.w[n] = cval1(i,j-q);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcslpara4_count;++g)
    {
    i=p->gcslpara4[g][0];
    j=p->gcslpara4[g][1];
        
        if(p->gcslpara4[g][3]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			C1.n[n] = cval1(i+q+2,j);
			C1.s[n] = cval1(i+q,j);
			}
		
		++n;
		}
	}
}

void ghostcell::column2D_pt2(lexer* p, sliceint &cval2, cpt2D &C2)
{
	n=0;
    SLICELOOP2
	{	
	C2.p[n] = cval2(i,j);
	
	C2.n[n] = cval2(i+1,j);
	C2.s[n] = cval2(i-1,j);
	
	C2.w[n] = cval2(i,j+1);
	C2.e[n] = cval2(i,j-1);

	++n;
	}
	
	GGCSL2LOOP
    {
    i=p->gcbsl2[g][0];
    j=p->gcbsl2[g][1];

        if(p->gcbsl2[g][3]==1)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C2.n[n] = cval2(i-q,j);
			C2.s[n] = cval2(i-q-2,j);
			}

        ++n;
        }

        if(p->gcbsl2[g][3]==2)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C2.e[n] = cval2(i,j+q);
			C2.w[n] = cval2(i,j+q+2);
			}
			
        ++n;
        }

        if(p->gcbsl2[g][3]==3)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C2.e[n] = cval2(i,j-q-2); 
			C2.w[n] = cval2(i,j-q); 
			}

        ++n;
        }

        if(p->gcbsl2[g][3]==4)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C2.n[n] = cval2(i+q+2,j);
			C2.s[n] = cval2(i+q,j);
			}

		++n;
        }
    }
	
	
	for(g=0;g<p->gcslpara1_count;++g)
    {
    i=p->gcslpara1[g][0];
    j=p->gcslpara1[g][1];
        
        if(p->gcslpara1[g][4]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			C2.n[n] = cval2(i-q,j);
			C2.s[n] = cval2(i-q-2,j);
			}
 
		++n;
		}
	}
	
	for(g=0;g<p->gcslpara2_count;++g)
    {
    i=p->gcslpara2[g][0];
    j=p->gcslpara2[g][1];
        
        if(p->gcslpara2[g][4]==1)
		for(q=0;q<margin;++q)   
		{
			if(q<margin-1)
			{
			C2.e[n] = cval2(i,j+q);
			C2.w[n] = cval2(i,j+q+2);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcslpara3_count;++g)
    {
    i=p->gcslpara3[g][0];
    j=p->gcslpara3[g][1];
        
        if(p->gcslpara3[g][4]==1)
		for(q=0;q<margin;++q)   
		{
			if(q<margin-1)
			{
			C2.e[n] = cval2(i,j-q-2);
			C2.w[n] = cval2(i,j-q);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcslpara4_count;++g)
    {
    i=p->gcslpara4[g][0];
    j=p->gcslpara4[g][1];
        
        if(p->gcslpara4[g][4]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			C2.n[n] = cval2(i+q+2,j);
			C2.s[n] = cval2(i+q,j);
			}
		
		++n;
		}
	}
}

void ghostcell::column2D_pt4(lexer* p, sliceint &cval4, cpt2D &C4)
{
	n=0;
    SLICELOOP4
	{	
	C4.p[n] = cval4(i,j);
	
	C4.n[n] = cval4(i+1,j);
	C4.s[n] = cval4(i-1,j);
	
	C4.w[n] = cval4(i,j+1);
	C4.e[n] = cval4(i,j-1);

	++n;
	}
	
	GGCSL4LOOP
    {
    i=p->gcbsl4[g][0];
    j=p->gcbsl4[g][1];

        if(p->gcbsl4[g][3]==1)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C4.n[n] = cval4(i-q,j);
			C4.s[n] = cval4(i-q-2,j);
			}

        ++n;
        }

        if(p->gcbsl4[g][3]==2)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C4.e[n] = cval4(i,j+q);
			C4.w[n] = cval4(i,j+q+2);
			}
			
        ++n;
        }

        if(p->gcbsl4[g][3]==3)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C4.e[n] = cval4(i,j-q-2); 
			C4.w[n] = cval4(i,j-q); 
			}

        ++n;
        }

        if(p->gcbsl4[g][3]==4)
        for(q=0;q<margin;++q)
        {
			if(q<margin-1)
			{
			C4.n[n] = cval4(i+q+2,j);
			C4.s[n] = cval4(i+q,j);
			}

		++n;
        }
    }
	
	for(g=0;g<p->gcslpara1_count;++g)
    {
    i=p->gcslpara1[g][0];
    j=p->gcslpara1[g][1];
        
        if(p->gcslpara1[g][6]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			C4.n[n] = cval4(i-q,j);
			C4.s[n] = cval4(i-q-2,j);
			}
 
		++n;
		}
	}
	
	for(g=0;g<p->gcslpara2_count;++g)
    {
    i=p->gcslpara2[g][0];
    j=p->gcslpara2[g][1];
        
        if(p->gcslpara2[g][6]==1)
		for(q=0;q<margin;++q)   
		{
			if(q<margin-1)
			{
			C4.e[n] = cval4(i,j+q);
			C4.w[n] = cval4(i,j+q+2);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcslpara3_count;++g)
    {
    i=p->gcslpara3[g][0];
    j=p->gcslpara3[g][1];
        
        if(p->gcslpara3[g][6]==1)
		for(q=0;q<margin;++q)   
		{
			if(q<margin-1)
			{
			C4.e[n] = cval4(i,j-q-2);
			C4.w[n] = cval4(i,j-q);
			}
		
		++n;
		}
	}
	
	for(g=0;g<p->gcslpara4_count;++g)
    {
    i=p->gcslpara4[g][0];
    j=p->gcslpara4[g][1];
        
        if(p->gcslpara4[g][6]==1)
		for(q=0;q<margin;++q)
		{
			if(q<margin-1)
			{
			C4.n[n] = cval4(i+q+2,j);
			C4.s[n] = cval4(i+q,j);
			}
		
		++n;
		}
	}
}