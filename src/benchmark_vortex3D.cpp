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

#include"benchmark_vortex3D.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"

benchmark_vortex3D::benchmark_vortex3D(lexer *p, fdm *a)
{
    double r,xc,yc,zc,radius;
    double H;

    xc = 0.35;
    yc = 0.35;
	zc = 0.35;
    radius = 0.15;

    LOOP
    a->phi(i,j,k)=-1.0;

	LOOP
	{
    r = sqrt( pow(p->pos_x()-xc,2.0) + pow(p->pos_y()-yc,2.0) + pow(p->pos_z()-zc,2.0));
	if(r<=radius)
	a->phi(i,j,k)=1.0;
	}

	
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

benchmark_vortex3D::~benchmark_vortex3D()
{
}

void benchmark_vortex3D::start(lexer* p, fdm *a, ghostcell *pgc, convection *pconvec )
{
    double xc,yc,zc;

    ULOOP
    {
    xc = p->pos_x() + 0.5*p->DXM;
    yc = p->pos_y();
	zc = p->pos_z();

    a->u(i,j,k) = 2.0* pow(sin(PI*xc),2.0) * sin(2.0*PI*yc) * sin(2.0*PI*zc) * cos((PI*p->simtime)/3.0);
    }

    VLOOP
    {
    xc = p->pos_x();
    yc = p->pos_y() + 0.5*p->DXM;
	zc = p->pos_z();

    a->v(i,j,k) = -sin(2.0*PI*xc) * pow(sin(PI*yc),2.0) * sin(2.0*PI*zc) * cos((PI*p->simtime)/3.0);
    }
	
	WLOOP
    {
    xc = p->pos_x();
    yc = p->pos_y();
	zc = p->pos_z() + 0.5*p->DXM;

    a->w(i,j,k) = -sin(2.0*PI*xc) * sin(2.0*PI*yc) * pow(sin(PI*zc),2.0) * cos((PI*p->simtime)/3.0);
    }

    pgc->start1(p,a->u,10);
    pgc->start2(p,a->v,11);
	pgc->start2(p,a->v,12);
}
