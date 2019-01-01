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
#include"fdm.h"
#include"cart1.h"
#include"cart2.h"
#include"cart3.h"
#include"cart4.h"
#include"sliceint.h"

int ghostcell::column2D_pt1_count(lexer* p, fdm2D *b)
{
    count=0;

    SLICELOOP1
    ++count;

	GC1LOOP
    {
    i=p->gcbsl1[n][0];
    j=p->gcbsl1[n][1];

        if(p->gcbsl1[n][3]==1)
        for(q=0;q<margin;++q)
        ++count;

        if(p->gcbsl1[n][3]==2)
        for(q=0;q<margin;++q)
        ++count;
        
        if(p->gcbsl1[n][3]==3)
        for(q=0;q<margin;++q)
        ++count;

        if(p->gcbsl1[n][3]==4)
        for(q=0;q<margin;++q)
        ++count;
    }
	
	for(n=0;n<p->gcslpara1_count;++n)
    {
    i=p->gcslpara1[n][0];
    j=p->gcslpara1[n][1];
        
        if(p->gcslpara1[n][6]==1)
		for(q=0;q<margin;++q)
        ++count;
	}
	
	for(n=0;n<p->gcslpara2_count;++n)
    {
    i=p->gcslpara2[n][0];
    j=p->gcslpara2[n][1];
        
        if(p->gcslpara2[n][6]==1)
		for(q=0;q<margin;++q)
        ++count;
	}
	
	for(n=0;n<p->gcslpara3_count;++n)
    {
    i=p->gcslpara3[n][0];
    j=p->gcslpara3[n][1];
        
        if(p->gcslpara3[n][6]==1)
		for(q=0;q<margin;++q)
        ++count;
	}
	
	for(n=0;n<p->gcslpara4_count;++n)
    {
    i=p->gcslpara4[n][0];
    j=p->gcslpara4[n][1];
    
        if(p->gcslpara4[n][6]==1)
		for(q=0;q<margin;++q)
        ++count;
	}
    
    return count;
}

int ghostcell::column2D_pt2_count(lexer* p, fdm2D *b)
{
    count=0;

    SLICELOOP2
    ++count;

	GC2LOOP
    {
    i=p->gcbsl2[n][0];
    j=p->gcbsl2[n][1];

        if(p->gcbsl2[n][3]==1)
        for(q=0;q<margin;++q)
        ++count;

        if(p->gcbsl2[n][3]==2)
        for(q=0;q<margin;++q)
        ++count;
        
        if(p->gcbsl2[n][3]==3)
        for(q=0;q<margin;++q)
        ++count;

        if(p->gcbsl2[n][3]==4)
        for(q=0;q<margin;++q)
        ++count;
    }
	
	for(n=0;n<p->gcslpara1_count;++n)
    {
    i=p->gcslpara1[n][0];
    j=p->gcslpara1[n][1];
        
        if(p->gcslpara1[n][6]==1)
		for(q=0;q<margin;++q)
        ++count;
	}
	
	for(n=0;n<p->gcslpara2_count;++n)
    {
    i=p->gcslpara2[n][0];
    j=p->gcslpara2[n][1];
        
        if(p->gcslpara2[n][6]==1)
		for(q=0;q<margin;++q)
        ++count;
	}
	
	for(n=0;n<p->gcslpara3_count;++n)
    {
    i=p->gcslpara3[n][0];
    j=p->gcslpara3[n][1];
        
        if(p->gcslpara3[n][6]==1)
		for(q=0;q<margin;++q)
        ++count;
	}
	
	for(n=0;n<p->gcslpara4_count;++n)
    {
    i=p->gcslpara4[n][0];
    j=p->gcslpara4[n][1];
    
        if(p->gcslpara4[n][6]==1)
		for(q=0;q<margin;++q)
        ++count;
	}
    
    return count;
}


int ghostcell::column2D_pt4_count(lexer* p, fdm2D *b)
{
    count=0;

    SLICELOOP4
    ++count;

	GC4LOOP
    {
    i=p->gcbsl4[n][0];
    j=p->gcbsl4[n][1];

        if(p->gcbsl4[n][3]==1)
        for(q=0;q<margin;++q)
        ++count;

        if(p->gcbsl4[n][3]==2)
        for(q=0;q<margin;++q)
        ++count;
        
        if(p->gcbsl4[n][3]==3)
        for(q=0;q<margin;++q)
        ++count;

        if(p->gcbsl4[n][3]==4)
        for(q=0;q<margin;++q)
        ++count;
    }
	
	for(n=0;n<p->gcslpara1_count;++n)
    {
    i=p->gcslpara1[n][0];
    j=p->gcslpara1[n][1];
        
        if(p->gcslpara1[n][6]==1)
		for(q=0;q<margin;++q)
        ++count;
	}
	
	for(n=0;n<p->gcslpara2_count;++n)
    {
    i=p->gcslpara2[n][0];
    j=p->gcslpara2[n][1];
        
        if(p->gcslpara2[n][6]==1)
		for(q=0;q<margin;++q)
        ++count;
	}
	
	for(n=0;n<p->gcslpara3_count;++n)
    {
    i=p->gcslpara3[n][0];
    j=p->gcslpara3[n][1];
        
        if(p->gcslpara3[n][6]==1)
		for(q=0;q<margin;++q)
        ++count;
	}
	
	for(n=0;n<p->gcslpara4_count;++n)
    {
    i=p->gcslpara4[n][0];
    j=p->gcslpara4[n][1];
    
        if(p->gcslpara4[n][6]==1)
		for(q=0;q<margin;++q)
        ++count;
	}
    
    return count;
}