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

void grid::cval_gcb1(lexer* p, fieldint &cval)
{
	GC1LOOP
    {
    i=p->gcb1[n][0];
    j=p->gcb1[n][1];
    k=p->gcb1[n][2];
	
	p->gcb1[n][5]=cval(i,j,k);
	}
}

void grid::cval_gcb2(lexer* p, fieldint &cval)
{
	GC2LOOP
    {
    i=p->gcb2[n][0];
    j=p->gcb2[n][1];
    k=p->gcb2[n][2];
	
	p->gcb2[n][5]=cval(i,j,k);
	}
}

void grid::cval_gcb3(lexer* p, fieldint &cval)
{
	GC3LOOP
    {
    i=p->gcb3[n][0];
    j=p->gcb3[n][1];
    k=p->gcb3[n][2];
	
	p->gcb3[n][5]=cval(i,j,k);
	}
}

void grid::cval_gcb4(lexer* p, fieldint &cval)
{
	GC4LOOP
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];
	
	p->gcb4[n][5]=cval(i,j,k);
	}
}

void grid::cval_gcb4a(lexer* p, fieldint &cval)
{
	GC4ALOOP
    {
    i=p->gcb4a[n][0];
    j=p->gcb4a[n][1];
    k=p->gcb4a[n][2];
	
	p->gcb4a[n][5]=cval(i,j,k);
	}
}

void grid::cval_gcb6(lexer* p, fieldint &cval)
{
    int count1,count2;
    count1=count2=0;
    
	GC6LOOP
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];
	
	p->gcb6[n]=cval(i,j,k);
    
        if(p->gcb4[n][4]==1 || p->gcb4[n][4]==6)
        {
        p->gcin6[count1][3] = p->gcb6[n];
        ++count1;
        }
        
        if(p->gcb4[n][4]==2 || p->gcb4[n][4]==7 || p->gcb4[n][4]==8)
        {
        p->gcout6[count2][3] = p->gcb6[n];
        ++count2;
        }
	}
}

