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
#include"initialize.h"
#include"fdm.h"
#include"lexer.h"
#include"field4.h"
#include"ghostcell.h"
#include"reinifluid_RK3.h"
#include"reini_RK3.h"


void initialize::inivofPLIC(fdm*a, lexer* p, ghostcell* pgc)
{

double dx=p->DXM;
double r;
double vofdiff, xdiff;
p->phimean=p->F56;


    LOOP
	a->vof(i,j,k)=0.0;


	LOOP
	if(double(i)*dx+p->originx>=p->F51 && double(i)*dx+p->originx<p->F54
	&& double(j)*dx+p->originy>=p->F52 && double(j)*dx+p->originy<p->F55
	&& double(k)*dx+p->originz>=p->F53 && double(k)*dx+p->originz<p->F56)
	a->vof(i,j,k)=1.0;
    
    LOOP
    if(a->vof(i,j,k)<0.001 && (a->vof(i-1,j,k)>0.999 || a->vof(i,j,k-1)>0.999))
    a->vof(i,j,k)=0.3;
    


if(p->F57_1>0||p->F57_2>0||p->F57_3>0||p->F57_4>0)
{
	LOOP
	if(p->F57_1*((double(i)+0.5)*dx + p->originx )+ p->F57_2*((double(j)+0.5)*dx + p->originy )+ p->F57_3*((double(k)+0.5)*dx + p->originz ) < p->F57_4)
	a->vof(i,j,k)=1.0;
}

if(p->F58_4>0.0)
{
    p->F58_1 -= p->originx;
    p->F58_2 -= p->originy;
    p->F58_3 -= p->originz;

	LOOP
	{
    r = sqrt( pow((double(i)+0.5)*dx-p->F58_1,2.0)+pow((double(j)+0.5)*dx-p->F58_2,2.0)+pow((double(k)+0.5)*dx-p->F58_3,2.0));
	if(r<=p->F58_4)
	a->vof(i,j,k)=1.0;
	}
}

if(p->F60>-1.0e20)
{
    LOOP
    a->vof(i,j,k)=p->F60-p->pos_z();

p->phimean=p->F60;
}

    if((p->F60>-1.0e20 || p->F56>-1.0e20) && p->F62>-1.0e-20&& p->F63>-1.0e-20  )
    {
        vofdiff=p->F62-p->phimean;
        xdiff=p->xcoormax-p->F63;

        LOOP
        if(p->pos_x() > p->F63)
        a->vof(i,j,k)=(vofdiff/xdiff)*(p->pos_x()-p->F63) + p->phimean    - p->pos_z() ;
    }

	double H=0.0;
        
    LOOP
        a->phi(i,j,k)=-0.5+a->vof(i,j,k);
        
    //preini->start(p,1)
    
    
    pgc->start4(p,a->phi,50);
    LOOP
	{
		if(a->phi(i,j,k)>(p->psi))
            H=1.0;

		if(a->phi(i,j,k)<-(p->psi))
            H=0.0;

		if(fabs(a->phi(i,j,k))<=(p->psi))
            H=0.5*(1.0 + a->phi(i,j,k)/(p->psi) + (1.0/PI)*sin((PI*a->phi(i,j,k))/(p->psi)));

		a->ro(i,j,k)= 0.0; p->W1*H + p->W3*(1.0-H);
		a->visc(i,j,k)= p->W2*H + p->W4*(1.0-H);
	}

   // redistance(a, p, pdisc, pgc, pflow, 20);
    
    

	pgc->start4(p,a->vof,50);
    
    
    pgc->start4(p,a->phi,50);
	pgc->start4(p,a->ro,1);
	pgc->start4(p,a->visc,1);
}

