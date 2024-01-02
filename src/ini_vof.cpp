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
#include"ghostcell.h"

void initialize::inivof(fdm*a, lexer* p, ghostcell* pgc)
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
	{
		H=a->vof(i,j,k);

		H=MAX(H,0.0);
		H=MIN(H,1.0);

		a->ro(i,j,k)= p->W1*H + p->W3*(1.0-H);
		a->visc(i,j,k)= p->W2*H + p->W4*(1.0-H);
	}

	pgc->start4(p,a->vof,50);
	pgc->start4(p,a->ro,1);
	pgc->start4(p,a->visc,1);
}

void initialize::inivof_io(fdm*a, lexer* p, ghostcell* pgc)
{

    if(p->F61>-1.0e20)
    GC4LOOP
    {
        if(p->gcb4[n][4]==1)
        {
        i=p->gcb4[n][0];
        j=p->gcb4[n][1];
        k=p->gcb4[n][2];

        a->vof(i-1,j,k)=p->F61-p->pos_z();
        a->vof(i-2,j,k)=p->F61-p->pos_z();
        a->vof(i-3,j,k)=p->F61-p->pos_z();
        }
    }

    if(p->F62>-1.0e20)
    GC4LOOP
    {
        if(p->gcb4[n][4]==2)
        {
        i=p->gcb4[n][0];
        j=p->gcb4[n][1];
        k=p->gcb4[n][2];

        a->vof(i+1,j,k)=p->F62-p->pos_z();
        a->vof(i+2,j,k)=p->F62-p->pos_z();
        a->vof(i+3,j,k)=p->F62-p->pos_z();
        }
    }
}

void initialize::inivof_box(lexer* p, fdm *a, ghostcell* pgc)
{
    int istart, iend, jstart, jend, kstart, kend;
    int qn;


    for(qn=0;qn<p->F70;++qn)
    {
        istart = conv((p->F70_xs[qn]-p->originx)/p->DXM);
        iend = conv((p->F70_xe[qn]-p->originx)/p->DXM);

        jstart = conv((p->F70_ys[qn]-p->originy)/p->DXM);
        jend = conv((p->F70_ye[qn]-p->originy)/p->DXM);

        kstart = conv((p->F70_zs[qn]-p->originz)/p->DXM);
        kend = conv((p->F70_ze[qn]-p->originz)/p->DXM);

        //cout<<p->mpirank<<" F70: "<<p->F70<<" . "<<istart<<" "<<iend<<" "<<jstart<<" "<<jend<<" "<<kstart<<" "<<kend<<endl;

        LOOP
        if(i>=istart && i<iend && j>=jstart && j<jend && k>=kstart && k<kend)
        a->vof(i,j,k)=0.0;
    }

}

