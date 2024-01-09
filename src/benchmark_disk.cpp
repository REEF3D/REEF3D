/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"benchmark_disk.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"

benchmark_disk::benchmark_disk(lexer *p, fdm *a)
{
    double r,xc,zc,radius;
    double H;
	double xs,xe,zs,ze;
	
    xc = 0.0;
    zc = 0.5;
    radius = 0.3;

	xs = -0.05;
	xe = 0.05;
	zs = 0.2;
	ze = 0.7;

    LOOP
	{
		a->vof(i,j,k) = 0.0;
		
		r = sqrt(pow(p->pos_x() - xc, 2.0) + pow(p->pos_z() - zc, 2.0));
		if (r <= radius)
		{
			a->vof(i,j,k) = 1.0;
		}

		if (p->pos_x() >= xs && p->pos_x() < xe && p->pos_z() >= zs && p->pos_z() < ze)
		{
			a->vof(i,j,k) = 0.0;
		}
		
		a->test(i,j,k) = a->vof(i,j,k);		
	}



	// Inverse field
	if(p->F151==1)
	LOOP
    a->phi(i,j,k)*=-1.0;

    LOOP
	{
		if(a->phi(i,j,k)>=p->F45*p->DXM)
		H=1.0;

		if(a->phi(i,j,k)<-p->F45*p->DXM)
		H=0.0;

		if(fabs(a->phi(i,j,k))<=p->F45*p->DXM)
		H=0.5*(1.0 + a->phi(i,j,k)/p->F45*p->DXM + (1.0/PI)*sin((PI*a->phi(i,j,k))/p->F45*p->DXM));

		a->ro(i,j,k)= p->W1*H + p->W3*(1.0-H);
		a->visc(i,j,k)= p->W2*H + p->W4*(1.0-H);
	}    
}

benchmark_disk::~benchmark_disk()
{
}

void benchmark_disk::start(lexer* p, fdm *a, ghostcell *pgc, convection *pconvec )
{
    LOOP
    {
		a->u(i,j,k) = -2.0*PI*p->pos_z();
		a->v(i,j,k) = 0.0;
		a->w(i,j,k) = 2.0*PI*p->pos_x();
    }

    pgc->start1(p,a->u,10);
    pgc->start2(p,a->v,11);
	pgc->start2(p,a->w,12);
}
