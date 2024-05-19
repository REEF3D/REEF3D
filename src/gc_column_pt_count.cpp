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

int ghostcell::column_pt4_count(lexer* p, fdm* a)
{
	n=0;
    LOOP
	++n;

	GGC4LOOP
    {
    i=p->gcb4[g][0];
    j=p->gcb4[g][1];
    k=p->gcb4[g][2];

        for(q=0;q<margin;++q)
        ++n;
    }
	
	for(g=0;g<p->gcpara1_count;++g)
    {
    i=p->gcpara1[g][0];
    j=p->gcpara1[g][1];
    k=p->gcpara1[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara1[g][3]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara2_count;++g)
    {
    i=p->gcpara2[g][0];
    j=p->gcpara2[g][1];
    k=p->gcpara2[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara2[g][3]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara3_count;++g)
    {
    i=p->gcpara3[g][0];
    j=p->gcpara3[g][1];
    k=p->gcpara3[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara3[g][3]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara4_count;++g)
    {
    i=p->gcpara4[g][0];
    j=p->gcpara4[g][1];
    k=p->gcpara4[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara4[g][3]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara5_count;++g)
    {
    i=p->gcpara5[g][0];
    j=p->gcpara5[g][1];
    k=p->gcpara5[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara5[g][3]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara6_count;++g)
    {
    i=p->gcpara6[g][0];
    j=p->gcpara6[g][1];
    k=p->gcpara6[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara6[g][3]==1)
		++n;
	}
    
    return n;
}

int ghostcell::column_pt4a_count(lexer* p, fdm* a)
{
	n=0;
    ALOOP
	++n;
	
	GGC4ALOOP
    {
    i=p->gcb4[g][0];
    j=p->gcb4[g][1];
    k=p->gcb4[g][2];

        for(q=0;q<margin;++q)
        ++n;
    }
	
	for(g=0;g<p->gcpara1_count;++g)
    {
    i=p->gcpara1[g][0];
    j=p->gcpara1[g][1];
    k=p->gcpara1[g][2];
	
		for(q=0;q<margin;++q)
		++n;
	}
	
	for(g=0;g<p->gcpara2_count;++g)
    {
    i=p->gcpara2[g][0];
    j=p->gcpara2[g][1];
    k=p->gcpara2[g][2];
	
		for(q=0;q<margin;++q)
		++n;
	}
	
	for(g=0;g<p->gcpara3_count;++g)
    {
    i=p->gcpara3[g][0];
    j=p->gcpara3[g][1];
    k=p->gcpara3[g][2];
	
		for(q=0;q<margin;++q)
		++n;
	}
	
	for(g=0;g<p->gcpara4_count;++g)
    {
    i=p->gcpara4[g][0];
    j=p->gcpara4[g][1];
    k=p->gcpara4[g][2];
	
		for(q=0;q<margin;++q)
		++n;
	}
	
	for(g=0;g<p->gcpara5_count;++g)
    {
    i=p->gcpara5[g][0];
    j=p->gcpara5[g][1];
    k=p->gcpara5[g][2];
	
		for(q=0;q<margin;++q)
		++n;
	}
	
	for(g=0;g<p->gcpara6_count;++g)
    {
    i=p->gcpara6[g][0];
    j=p->gcpara6[g][1];
    k=p->gcpara6[g][2];
	
		for(q=0;q<margin;++q)
		++n;
	}
    
    return n;
}

int ghostcell::column_pt6_count(lexer* p, fdm* a)
{
	n=0;
    BASELOOP
	++n;

	GGC6LOOP
    {
    i=p->gcb4[g][0];
    j=p->gcb4[g][1];
    k=p->gcb4[g][2];

        for(q=0;q<margin;++q)
        ++n;
    }

	for(g=0;g<p->gcpara1_count;++g)
    {
    i=p->gcpara1[g][0];
    j=p->gcpara1[g][1];
    k=p->gcpara1[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara1[g][8]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara2_count;++g)
    {
    i=p->gcpara2[g][0];
    j=p->gcpara2[g][1];
    k=p->gcpara2[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara2[g][8]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara3_count;++g)
    {
    i=p->gcpara3[g][0];
    j=p->gcpara3[g][1];
    k=p->gcpara3[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara3[g][8]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara4_count;++g)
    {
    i=p->gcpara4[g][0];
    j=p->gcpara4[g][1];
    k=p->gcpara4[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara4[g][8]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara5_count;++g)
    {
    i=p->gcpara5[g][0];
    j=p->gcpara5[g][1];
    k=p->gcpara5[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara5[g][8]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara6_count;++g)
    {
    i=p->gcpara6[g][0];
    j=p->gcpara6[g][1];
    k=p->gcpara6[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara6[g][8]==1)
		++n;
	}
    
    return n;
}
