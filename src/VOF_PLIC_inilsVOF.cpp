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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"VOF_PLIC.h"
#include"fdm.h"
#include"lexer.h"
#include"ghostcell.h"

void VOF_PLIC::iniphi(fdm*a, lexer* p, ghostcell* pgc)
{

double dx=p->DXM;
double r;
double phidiff, xdiff;
p->phimean=p->F56;


    LOOP
	a->phi(i,j,k)=-1.0;
	
	pgc->start4(p,a->phi,50);

    if(p->F50_flag==1)
	LOOP
	if(double(i)*dx+p->originx>=p->F51 && double(i)*dx+p->originx<p->F54
	&& double(j)*dx+p->originy>=p->F52 && double(j)*dx+p->originy<p->F55
	&& double(k)*dx+p->originz>=p->F53 && double(k)*dx+p->originz<p->F56)
	a->phi(i,j,k)=1.0;


if(p->F57_1>0||p->F57_2>0||p->F57_3>0||p->F57_4>0)
{
	LOOP
	if(p->F57_1*((double(i)+0.5)*dx + p->originx )+ p->F57_2*((double(j)+0.5)*dx + p->originy )+ p->F57_3*((double(k)+0.5)*dx + p->originz ) < p->F57_4)
	a->phi(i,j,k)=1.0;
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
	a->phi(i,j,k)=1.0;
	}
}

if(p->F59_r>0.0)
{
    
    p->F59_xm -= p->originx;
    p->F59_ym -= p->originy;
	
	LOOP
	{
    r = sqrt( pow((double(i)+0.5)*dx-p->F59_xm,2.0)+pow((double(j)+0.5)*dx-p->F59_ym,2.0));

	if(r<=p->F59_r && p->pos_z()>p->F59_zs && p->pos_z()<=p->F59_ze)
	a->phi(i,j,k)=1.0;
	}
}

if(p->F60>-1.0e20)
{
    LOOP
    a->phi(i,j,k)=p->F60-p->pos_z();

	p->phimean=p->F60;
	
}
    double epsi = p->F45*p->DXM;
    

    if((p->F60>-1.0e20 || p->F56>-1.0e20) && p->F62>-1.0e-20&& p->F63>-1.0e-20  )
    {
        phidiff=p->F62-p->phimean;
        xdiff=p->xcoormax-p->F63;

        LOOP
        if(p->pos_x() > p->F63)
        a->phi(i,j,k)=(phidiff/xdiff)*(p->pos_x()-p->F63) + p->phimean    - p->pos_z() ;
    }

	double H=0.0;

	LOOP
	{
		if(a->phi(i,j,k)>=0)
		H=1.0;

		if(a->phi(i,j,k)<0)
		H=0.0;

		if(fabs(a->phi(i,j,k))<=epsi)
		H=0.5*(1.0 + a->phi(i,j,k)/epsi + (1.0/PI)*sin((PI*a->phi(i,j,k))/epsi));

		a->ro(i,j,k)= p->W1*H + p->W3*(1.0-H);
		a->visc(i,j,k)= p->W2*H + p->W4*(1.0-H);
	}

	pgc->start4(p,a->phi,50);
	pgc->start4(p,a->ro,1);
	pgc->start4(p,a->visc,1);


	p->phimean=p->F56;

    if(p->F60>-1.0e20)
    {
    p->phimean=p->F60;
    p->phiout=p->F60;
    }

    if(p->F62>-1.0e20)
    p->phiout=p->F62;

    if(p->F61>-1.0e20)
    p->phiin=p->F62;
	
	
	if(p->F64==1)
	LOOP
	{
	
	//a->phi(i,j,k) = p->F61-p->pos_z()
	}
	
	
	pgc->start4(p,a->phi,50);
}

void VOF_PLIC::iniphi_io(fdm*a, lexer* p, ghostcell* pgc)
{
    if(p->F61>-1.0e20)
    GC4LOOP
    {
        if(p->gcb4[n][4]==1)
        {
        i=p->gcb4[n][0];
        j=p->gcb4[n][1];
        k=p->gcb4[n][2];

        a->phi(i-1,j,k)=p->F61-p->pos_z();
        a->phi(i-2,j,k)=p->F61-p->pos_z();
        a->phi(i-3,j,k)=p->F61-p->pos_z();
        }
    p->phiin=p->F61;
    }

    if(p->F62>-1.0e20)
    GC4LOOP
    {
        if(p->gcb4[n][4]==2)
        {
        i=p->gcb4[n][0];
        j=p->gcb4[n][1];
        k=p->gcb4[n][2];

        a->phi(i+1,j,k)=p->F62-p->pos_z();
        a->phi(i+2,j,k)=p->F62-p->pos_z();
        a->phi(i+3,j,k)=p->F62-p->pos_z();
        }
    } 
}

void VOF_PLIC::iniphi_box(lexer* p, fdm *a, ghostcell* pgc)
{
    int istart, iend, jstart, jend, kstart, kend;
    int qn;
    
    if(p->F70>0)
    LOOP
	a->phi(i,j,k)=-1.0;

    for(qn=0;qn<p->F70;++qn)
    {
        istart = conv((p->F70_xs[qn]-p->originx)/p->DXM);
        iend = conv((p->F70_xe[qn]-p->originx)/p->DXM);

        jstart = conv((p->F70_ys[qn]-p->originy)/p->DXM);
        jend = conv((p->F70_ye[qn]-p->originy)/p->DXM);

        kstart = conv((p->F70_zs[qn]-p->originz)/p->DXM);
        kend = conv((p->F70_ze[qn]-p->originz)/p->DXM);


        LOOP
        if(i>=istart && i<iend && j>=jstart && j<jend && k>=kstart && k<kend)
        a->phi(i,j,k)=1;
    }
	
	for(qn=0;qn<p->F71;++qn)
    {
        istart = conv((p->F71_xs[qn]-p->originx)/p->DXM);
        iend = conv((p->F71_xe[qn]-p->originx)/p->DXM);

        jstart = conv((p->F71_ys[qn]-p->originy)/p->DXM);
        jend = conv((p->F71_ye[qn]-p->originy)/p->DXM);

        kstart = conv((p->F71_zs[qn]-p->originz)/p->DXM);
        kend = conv((p->F71_ze[qn]-p->originz)/p->DXM);


        LOOP
        if(i>=istart && i<iend && j>=jstart && j<jend && k>=kstart && k<kend)
        a->phi(i,j,k)=-1;
    }
	
	for(qn=0;qn<p->F72;++qn)
    {
		istart = conv((p->F72_xs[qn]-p->originx)/p->DXM);
        iend = conv((p->F72_xe[qn]-p->originx)/p->DXM);

        jstart = conv((p->F72_ys[qn]-p->originy)/p->DXM);
        jend = conv((p->F72_ye[qn]-p->originy)/p->DXM);
		

        LOOP
        if(i>=istart && i<iend && j>=jstart && j<jend)
        a->phi(i,j,k)=p->F72_h[qn]-p->pos_z();

	}
}

void VOF_PLIC::iniphi_surfarea(lexer* p, fdm *a, ghostcell* pgc)
{
	double dx,dy,dz,dnorm,dirac;
	double area=0.0;
	double epsi = 1.6*p->DXM;
	
    LOOP
	{
	dx = (a->phi(i+1,j,k)-a->phi(i-1,j,k))/(2.0*p->DXM);
	dy = (a->phi(i,j+1,k)-a->phi(i,j-1,k))/(2.0*p->DXM);
	dz = (a->phi(i,j,k+1)-a->phi(i,j,k-1))/(2.0*p->DXM);
	
	dnorm = sqrt(dx*dx + dy*dy + dz*dz);
	
	//if(fabs(a->phi(i,j,k))<epsi)
	//cout<<" dx: "<<dx<<" dy: "<<dy<<" dz: "<<dz<<" dnorm: "<<dnorm<<endl;
	
	dirac=0.0;
	
	if(fabs(a->phi(i,j,k))<epsi)
	dirac = (0.5/epsi)*(1.0 + cos((PI*a->phi(i,j,k))/epsi));
	
	area +=  pow(p->DXM,3.0) * dirac *dnorm;
	}
	
	area = pgc->globalsum(area);
	
	//if(p->mpirank==0)
	//cout<<"Surface Area: "<<area<<endl;

}


int VOF_PLIC::conv(double a)
{

	int b,c;
	double d,diff;

	c= int( a);
	d=double(c);
	diff=a-d;

	b=c;

	if(diff>0.5)
	b=c+1;

	if(diff<=-0.5)
	b=c-1;


	return b;

}
