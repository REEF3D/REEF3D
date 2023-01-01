/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"benchmark_convection.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"

benchmark_convection::benchmark_convection(lexer *p, fdm *a)
{
	LOOP
	a->phi(i,j,k)=0.0;

	LOOP
	{
		if(p->pos_x()>0.5 && p->pos_x()<1.5)
		a->phi(i,j,k) = 1.0;
		
		if(p->pos_x()>2.5 && p->pos_x()<3.5)
		a->phi(i,j,k) = 1.0 - fabs(2.0*(p->pos_x()-3.0));
		
		if(p->pos_x()>4.5 && p->pos_x()<5.5)
		a->phi(i,j,k) = sin((p->pos_x()-4.5)*PI);
	}
	
	
}

benchmark_convection::~benchmark_convection()
{
}

void benchmark_convection::start(lexer* p, fdm *a, ghostcell *pgc, convection *pconvec )
{
	if(p->count==0)
	LOOP
	{
		if(p->pos_x()>0.5 && p->pos_x()<1.5)
		a->phi(i,j,k) = 1.0;
		
		if(p->pos_x()>2.5 && p->pos_x()<3.5)
		a->phi(i,j,k) = 1.0 - fabs(2.0*(p->pos_x()-3.0));
		
		if(p->pos_x()>4.5 && p->pos_x()<5.5)
		a->phi(i,j,k) = sin((p->pos_x()-4.5)*PI);
	}
	
	 pgc->start4(p,a->phi,40);	
	 
    ULOOP
    a->u(i,j,k) = 1.0;

    pgc->start1(p,a->u,10);	
}
