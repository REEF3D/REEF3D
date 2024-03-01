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
#include"cart1.h"
#include"cart2.h"
#include"cart3.h"
#include"cart4.h"

void ghostcell::cval_gcpara4(lexer* p, fdm* a, fieldint &cval4)
{
	for(n=0;n<p->gcpara1_count;++n)
    if(p->gcpara1[n][6]==1)
    {
    i=p->gcpara1[n][0];
    j=p->gcpara1[n][1];
    k=p->gcpara1[n][2];
		
	p->gcpara1[n][12]=cval4(i,j,k);		
	}
	
	for(n=0;n<p->gcpara2_count;++n)
    if(p->gcpara2[n][6]==1)
    {
    i=p->gcpara2[n][0];
    j=p->gcpara2[n][1];
    k=p->gcpara2[n][2];
	
	p->gcpara2[n][12]=cval4(i,j,k);		
	}
	
	for(n=0;n<p->gcpara3_count;++n)
    if(p->gcpara3[n][6]==1)
    {
    i=p->gcpara3[n][0];
    j=p->gcpara3[n][1];
    k=p->gcpara3[n][2];
	
	p->gcpara3[n][12]=cval4(i,j,k);	
	}
	
	for(n=0;n<p->gcpara4_count;++n)
    if(p->gcpara4[n][6]==1)
    {
    i=p->gcpara4[n][0];
    j=p->gcpara4[n][1];
    k=p->gcpara4[n][2];
	
	p->gcpara4[n][12]=cval4(i,j,k);		
	}
	
	for(n=0;n<p->gcpara5_count;++n)
    if(p->gcpara5[n][6]==1)
    {
    i=p->gcpara5[n][0];
    j=p->gcpara5[n][1];
    k=p->gcpara5[n][2];

	p->gcpara5[n][12]=cval4(i,j,k);		
	}
	
	for(n=0;n<p->gcpara6_count;++n)
    if(p->gcpara6[n][6]==1)
    {
    i=p->gcpara6[n][0];
    j=p->gcpara6[n][1];
    k=p->gcpara6[n][2];

	p->gcpara6[n][12]=cval4(i,j,k);		
	}
}

void ghostcell::cval_gcpara4a(lexer* p, fdm* a, fieldint &cval4a)
{
	for(n=0;n<p->gcpara1_count;++n)
    if(p->gcpara1[n][7]==1)
    {
    i=p->gcpara1[n][0];
    j=p->gcpara1[n][1];
    k=p->gcpara1[n][2];
		
	p->gcpara1[n][13]=cval4a(i,j,k);		
	}
	
	for(n=0;n<p->gcpara2_count;++n)
    if(p->gcpara2[n][7]==1)
    {
    i=p->gcpara2[n][0];
    j=p->gcpara2[n][1];
    k=p->gcpara2[n][2];
	
	p->gcpara2[n][13]=cval4a(i,j,k);		
	}
	
	for(n=0;n<p->gcpara3_count;++n)
    if(p->gcpara3[n][7]==1)
    {
    i=p->gcpara3[n][0];
    j=p->gcpara3[n][1];
    k=p->gcpara3[n][2];
	
	p->gcpara3[n][13]=cval4a(i,j,k);		
	}
	
	for(n=0;n<p->gcpara4_count;++n)
    if(p->gcpara4[n][7]==1)
    {
    i=p->gcpara4[n][0];
    j=p->gcpara4[n][1];
    k=p->gcpara4[n][2];
	
	p->gcpara4[n][13]=cval4a(i,j,k);		
	}
	
	for(n=0;n<p->gcpara5_count;++n)
    if(p->gcpara5[n][7]==1)
    {
    i=p->gcpara5[n][0];
    j=p->gcpara5[n][1];
    k=p->gcpara5[n][2];

	p->gcpara5[n][13]=cval4a(i,j,k);		
	}
	
	for(n=0;n<p->gcpara6_count;++n)
    if(p->gcpara6[n][7]==1)
    {
    i=p->gcpara6[n][0];
    j=p->gcpara6[n][1];
    k=p->gcpara6[n][2];

	p->gcpara6[n][13]=cval4a(i,j,k);		
	}
}

void ghostcell::cval_gcpara6(lexer* p, fdm* a, fieldint &cval6)
{
	for(n=0;n<p->gcpara1_count;++n)
    {
    i=p->gcpara1[n][0];
    j=p->gcpara1[n][1];
    k=p->gcpara1[n][2];
		
	p->gcpara1[n][14]=cval6(i,j,k);		
	}
	
	for(n=0;n<p->gcpara2_count;++n)
    {
    i=p->gcpara2[n][0];
    j=p->gcpara2[n][1];
    k=p->gcpara2[n][2];
	
	p->gcpara2[n][14]=cval6(i,j,k);		
	}
	
	for(n=0;n<p->gcpara3_count;++n)
    {
    i=p->gcpara3[n][0];
    j=p->gcpara3[n][1];
    k=p->gcpara3[n][2];
	
	p->gcpara3[n][14]=cval6(i,j,k);		
	}
	
	for(n=0;n<p->gcpara4_count;++n)
    {
    i=p->gcpara4[n][0];
    j=p->gcpara4[n][1];
    k=p->gcpara4[n][2];
	
	p->gcpara4[n][14]=cval6(i,j,k);		
	}
	
	for(n=0;n<p->gcpara5_count;++n)
    {
    i=p->gcpara5[n][0];
    j=p->gcpara5[n][1];
    k=p->gcpara5[n][2];

	p->gcpara5[n][14]=cval6(i,j,k);		
	}
	
	for(n=0;n<p->gcpara6_count;++n)
    {
    i=p->gcpara6[n][0];
    j=p->gcpara6[n][1];
    k=p->gcpara6[n][2];

	p->gcpara6[n][14]=cval6(i,j,k);		
	}
}

