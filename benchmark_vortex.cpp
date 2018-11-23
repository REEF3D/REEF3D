/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
--------------------------------------------------------------------*/

#include"benchmark_vortex.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"

benchmark_vortex::benchmark_vortex(lexer *p, fdm *a)
{
    double xc,yc,radius;
	double dist,sign;
    double H;

    xc = 0.5;
    yc = 0.75;
    radius = 0.15;

    LOOP
    a->phi(i,j,k)=1.0;
/*
	LOOP
	{
    r = sqrt( pow(p->pos_x()-xc,2.0) + pow(p->pos_y()-yc,2.0));
	if(r<=radius)
	a->phi(i,j,k)=-1.0;
	}
	
	*/
	LOOP
	{
		
	dist = sqrt( pow(p->pos_x()-xc,2.0) + pow(p->pos_y()-yc,2.0));
	
	if(dist<=radius)
	sign=1.0;
	
	if(dist>radius)
	sign=-1.0;
	
	a->phi(i,j,k)=sign*dist;	
	}
	
	if(p->F151==1)
	LOOP
    a->phi(i,j,k)*=-1.0;

    LOOP
	{
		if(a->phi(i,j,k)>=p->F45*p->dx)
		H=1.0;

		if(a->phi(i,j,k)<-p->F45*p->dx)
		H=0.0;

		if(fabs(a->phi(i,j,k))<=p->F45*p->dx)
		H=0.5*(1.0 + a->phi(i,j,k)/p->F45*p->dx + (1.0/PI)*sin((PI*a->phi(i,j,k))/p->F45*p->dx));

		a->ro(i,j,k)= p->W1*H + p->W3*(1.0-H);
		a->visc(i,j,k)= p->W2*H + p->W4*(1.0-H);
	}
}

benchmark_vortex::~benchmark_vortex()
{
}

void benchmark_vortex::start(lexer* p, fdm *a, ghostcell *pgc, convection *pconvec )
{
    double xc,yc;

    ULOOP
    {
    xc = p->pos_x() + 0.5*p->dx;
    yc = p->pos_y();

    a->u(i,j,k) = -pow(sin(PI*xc),2.0) * sin(2.0*PI*yc) * cos((PI*p->simtime)/8.0);
    }

    VLOOP
    {
    xc = p->pos_x();
    yc = p->pos_y() + 0.5*p->dx;

    a->v(i,j,k) = pow(sin(PI*yc),2.0) * sin(2.0*PI*xc) * cos((PI*p->simtime)/8.0);
    }

    pgc->start1(p,a->u,10);
    pgc->start2(p,a->v,11);



}
