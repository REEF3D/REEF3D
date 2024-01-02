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

#include"sflow_filter.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"

sflow_filter::sflow_filter(lexer* p) : f1x(p), f1y(p), f2x(p), f2y(p), f4x(p), f4y(p) 
{
}

sflow_filter::~sflow_filter()
{
}

void sflow_filter::filter(lexer* p, fdm2D *b, ghostcell *pgc)
{
    if(p->count%200==0)
    {
    filter1(p,b,pgc);
    filter4(p,b,pgc);
    }
}

void sflow_filter::filter1(lexer* p, fdm2D *b, ghostcell *pgc)
{
	//SLICELOOP1
	//f1x(i,j) = (1.0/254.0)*(186.0*b->P(i,j) + 56.0*(b->P(i+1,j) + b->P(i-1,j)) - 28.0*(b->P(i+2,j) + b->P(i-2,j)) + 8.0*(b->P(i+3,j) + b->P(i-3,j)));
    
    SLICELOOP1
	f1x(i,j) = (1.0/4.0)*(2.0*b->P(i,j) + (b->P(i+1,j) + b->P(i-1,j)));
    
    //SLICELOOP1
	//f1y(i,j) = (1.0/254.0)*(186.0*f1x(i,j) + 56.0*(f1x(i,j+1) + f1x(i,j+1)) - 28.0*(f1x(i,j+2) + f1x(i,j-2)) + 8.0*(f1x(i,j+3) + f1x(i,j-3)));

    SLICELOOP1
    b->P(i,j) = f1x(i,j);
    
    pgc->gcsl_start1(p,b->P,10);
}

void sflow_filter::filter2(lexer* p, fdm2D *b, ghostcell *pgc)
{
}

void sflow_filter::filter4(lexer* p, fdm2D *b, ghostcell *pgc)
{
    //SLICELOOP4
	//f4x(i,j) = (1.0/254.0)*(186.0*b->eta(i,j) + 56.0*(b->eta(i+1,j) + b->eta(i-1,j)) - 28.0*(b->eta(i+2,j) + b->eta(i-2,j)) + 8.0*(b->eta(i+3,j) + b->eta(i-3,j)));
    
    SLICELOOP4
	f4x(i,j) = (1.0/4.0)*(2.0*b->eta(i,j) + (b->eta(i+1,j) + b->eta(i-1,j)));
    
    //SLICELOOP1
	//f1y(i,j) = (1.0/254.0)*(186.0*f1x(i,j) + 56.0*(f1x(i,j+1) + f1x(i,j+1)) - 28.0*(f1x(i,j+2) + f1x(i,j-2)) + 8.0*(f1x(i,j+3) + f1x(i,j-3)));

    SLICELOOP4
    b->eta(i,j) = f4x(i,j);
    
    pgc->gcsl_start4(p,b->eta,50);

}

