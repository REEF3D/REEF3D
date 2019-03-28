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
#include"fdm2D.h"

void ghostcell::cval_gcslpara1(lexer* p, fdm2D* b, sliceint &cval1)
{
	for(n=0;n<p->gcslpara1_count;++n)
    if(p->gcslpara1[n][3]==1)
    {
    i=p->gcslpara1[n][0];
    j=p->gcslpara1[n][1];
		
	p->gcslpara1[n][9]=cval1(i,j);		
	}
	
	for(n=0;n<p->gcslpara2_count;++n)
    if(p->gcslpara2[n][3]==1)
    {
    i=p->gcslpara2[n][0];
    j=p->gcslpara2[n][1];
	
	p->gcslpara2[n][9]=cval1(i,j);		
	}
	
	for(n=0;n<p->gcslpara3_count;++n)
    if(p->gcslpara3[n][3]==1)
    {
    i=p->gcslpara3[n][0];
    j=p->gcslpara3[n][1];
	
	p->gcslpara3[n][9]=cval1(i,j);		
	}
	
	for(n=0;n<p->gcslpara4_count;++n)
    if(p->gcslpara4[n][3]==1)
    {
    i=p->gcslpara4[n][0];
    j=p->gcslpara4[n][1];

	p->gcslpara4[n][9]=cval1(i,j);		
	}
}

void ghostcell::cval_gcslpara2(lexer* p, fdm2D* b, sliceint &cval2)
{
	for(n=0;n<p->gcslpara1_count;++n)
    if(p->gcslpara1[n][4]==1)
    {
    i=p->gcslpara1[n][0];
    j=p->gcslpara1[n][1];
		
	p->gcslpara1[n][10]=cval2(i,j);		
	}
	
	for(n=0;n<p->gcslpara2_count;++n)
    if(p->gcslpara2[n][4]==1)
    {
    i=p->gcslpara2[n][0];
    j=p->gcslpara2[n][1];
	
	p->gcslpara2[n][10]=cval2(i,j);		
	}
	
	for(n=0;n<p->gcslpara3_count;++n)
    if(p->gcslpara3[n][4]==1)
    {
    i=p->gcslpara3[n][0];
    j=p->gcslpara3[n][1];
	
	p->gcslpara3[n][10]=cval2(i,j);		
	}
	
	for(n=0;n<p->gcslpara4_count;++n)
    if(p->gcslpara4[n][4]==1)
    {
    i=p->gcslpara4[n][0];
    j=p->gcslpara4[n][1];
	
	p->gcslpara4[n][10]=cval2(i,j);		
	}
}

void ghostcell::cval_gcslpara4(lexer* p, fdm2D* b, sliceint &cval4)
{
	for(n=0;n<p->gcslpara1_count;++n)
    if(p->gcslpara1[n][6]==1)
    {
    i=p->gcslpara1[n][0];
    j=p->gcslpara1[n][1];
		
	p->gcslpara1[n][11]=cval4(i,j);		
	}
	
	for(n=0;n<p->gcslpara2_count;++n)
    if(p->gcslpara2[n][6]==1)
    {
    i=p->gcslpara2[n][0];
    j=p->gcslpara2[n][1];
	
	p->gcslpara2[n][11]=cval4(i,j);		
	}
	
	for(n=0;n<p->gcslpara3_count;++n)
    if(p->gcslpara3[n][6]==1)
    {
    i=p->gcslpara3[n][0];
    j=p->gcslpara3[n][1];
	
	p->gcslpara3[n][11]=cval4(i,j);	
	}
	
	for(n=0;n<p->gcslpara4_count;++n)
    if(p->gcslpara4[n][6]==1)
    {
    i=p->gcslpara4[n][0];
    j=p->gcslpara4[n][1];
	
	p->gcslpara4[n][11]=cval4(i,j);		
	}
	

}
