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
along with this program; if not, see <http:www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"
#include"cart1.h"
#include"cart2.h"
#include"cart3.h"
#include"cart4.h"
#include"sliceint.h"

void ghostcell::cval2Dupdate1(lexer* p, sliceint &cval1)
{
    count=0;

    SLICELOOP1
	{
    cval1(i,j)=count;
    ++count;
	}
	
	GCSL1LOOP
    {
    i=p->gcbsl1[n][0];
    j=p->gcbsl1[n][1];

        if(p->gcbsl1[n][3]==1)
        for(q=0;q<margin;++q)
        {
        cval1(i-1-q,j)=count;
        ++count;
        }

        if(p->gcbsl1[n][3]==2)
        for(q=0;q<margin;++q)
        {
        cval1(i,j+1+q)=count;
        ++count;
        }

        if(p->gcbsl1[n][3]==3)
        for(q=0;q<margin;++q)
        {
        cval1(i,j-1-q)=count;
        ++count;
        }

        if(p->gcbsl1[n][3]==4)
        for(q=0;q<margin;++q)
        {
        cval1(i+1+q,j)=count;
        ++count;
        }
    }
	
	for(n=0;n<p->gcslpara1_count;++n)
    {
    i=p->gcslpara1[n][0];
    j=p->gcslpara1[n][1];
        
        if(p->gcslpara1[n][3]==1)
		for(q=0;q<margin;++q)
        {
        cval1(i-1-q,j)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcslpara2_count;++n)
    {
    i=p->gcslpara2[n][0];
    j=p->gcslpara2[n][1];
        
        if(p->gcslpara2[n][3]==1)
		for(q=0;q<margin;++q)
        {
        cval1(i,j+1+q)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcslpara3_count;++n)
    {
    i=p->gcslpara3[n][0];
    j=p->gcslpara3[n][1];
        
        if(p->gcslpara3[n][3]==1)
		for(q=0;q<margin;++q)
        {
        cval1(i,j-1-q)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcslpara4_count;++n)
    {
    i=p->gcslpara4[n][0];
    j=p->gcslpara4[n][1];
    
        if(p->gcslpara4[n][3]==1)
		for(q=0;q<margin;++q)
        {
        cval1(i+1+q,j)=count;
        ++count;
        }
	}
}


void ghostcell::cval2Dupdate2(lexer* p, sliceint &cval2)
{
    count=0;

    SLICELOOP2
	{
    cval2(i,j)=count;
    ++count;
	}
	

	GCSL2LOOP
    {
    i=p->gcbsl2[n][0];
    j=p->gcbsl2[n][1];

        if(p->gcbsl2[n][3]==1)
        for(q=0;q<margin;++q)
        {
        cval2(i-1-q,j)=count;
        ++count;
        }

        if(p->gcbsl2[n][3]==2)
        for(q=0;q<margin;++q)
        {
        cval2(i,j+1+q)=count;
        ++count;
        }

        if(p->gcbsl2[n][3]==3)
        for(q=0;q<margin;++q)
        {
        cval2(i,j-1-q)=count;
        ++count;
        }

        if(p->gcbsl2[n][3]==4)
        for(q=0;q<margin;++q)
        {
        cval2(i+1+q,j)=count;
        ++count;
        }
    }
	
	for(n=0;n<p->gcslpara1_count;++n)
    {
    i=p->gcslpara1[n][0];
    j=p->gcslpara1[n][1];
        
        if(p->gcslpara1[n][4]==1)
		for(q=0;q<margin;++q)
        {
        cval2(i-1-q,j)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcslpara2_count;++n)
    {
    i=p->gcslpara2[n][0];
    j=p->gcslpara2[n][1];
        
        if(p->gcslpara2[n][4]==1)
		for(q=0;q<margin;++q)
        {
        cval2(i,j+1+q)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcslpara3_count;++n)
    {
    i=p->gcslpara3[n][0];
    j=p->gcslpara3[n][1];
        
        if(p->gcslpara3[n][4]==1)
		for(q=0;q<margin;++q)
        {
        cval2(i,j-1-q)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcslpara4_count;++n)
    {
    i=p->gcslpara4[n][0];
    j=p->gcslpara4[n][1];
    
        if(p->gcslpara4[n][4]==1)
		for(q=0;q<margin;++q)
        {
        cval2(i+1+q,j)=count;
        ++count;
        }
	}
}

void ghostcell::cval2Dupdate4(lexer* p, sliceint &cval4)
{
    count=0;

    SLICELOOP4
	{
    cval4(i,j)=count;
    ++count;
	}
	

	GCSL4LOOP
    {
    i=p->gcbsl4[n][0];
    j=p->gcbsl4[n][1];

        if(p->gcbsl4[n][3]==1)
        for(q=0;q<margin;++q)
        {
        cval4(i-1-q,j)=count;
        ++count;
        }

        if(p->gcbsl4[n][3]==2)
        for(q=0;q<margin;++q)
        {
        cval4(i,j+1+q)=count;
        ++count;
        }

        if(p->gcbsl4[n][3]==3)
        for(q=0;q<margin;++q)
        {
        cval4(i,j-1-q)=count;
        ++count;
        }

        if(p->gcbsl4[n][3]==4)
        for(q=0;q<margin;++q)
        {
        cval4(i+1+q,j)=count;
        ++count;
        }
    }
	
	for(n=0;n<p->gcslpara1_count;++n)
    {
    i=p->gcslpara1[n][0];
    j=p->gcslpara1[n][1];
        
        if(p->gcslpara1[n][6]==1)
		for(q=0;q<margin;++q)
        {
        cval4(i-1-q,j)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcslpara2_count;++n)
    {
    i=p->gcslpara2[n][0];
    j=p->gcslpara2[n][1];
        
        if(p->gcslpara2[n][6]==1)
		for(q=0;q<margin;++q)
        {
        cval4(i,j+1+q)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcslpara3_count;++n)
    {
    i=p->gcslpara3[n][0];
    j=p->gcslpara3[n][1];
        
        if(p->gcslpara3[n][6]==1)
		for(q=0;q<margin;++q)
        {
        cval4(i,j-1-q)=count;
        ++count;
        }
	}
	
	for(n=0;n<p->gcslpara4_count;++n)
    {
    i=p->gcslpara4[n][0];
    j=p->gcslpara4[n][1];
    
        if(p->gcslpara4[n][6]==1)
		for(q=0;q<margin;++q)
        {
        cval4(i+1+q,j)=count;
        ++count;
        }
	}
}