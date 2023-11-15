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
	if(p->XN[IP]>=p->F51 && p->XN[IP]<p->F54
	&& p->YN[JP]>=p->F52 && p->YN[JP]<p->F55
	&& p->ZN[KP]>=p->F53 && p->ZN[KP]<p->F56)
	a->phi(i,j,k)=1.0;


if(p->F57_1>0||p->F57_2>0||p->F57_3>0||p->F57_4>0)
{
	LOOP
	if(p->F57_1*p->XP[IP]+ p->F57_2*p->YP[JP]+ p->F57_3*p->ZP[KP] < p->F57_4)
	a->phi(i,j,k)=1.0;
}

if(p->F58_4>0.0)
{
	LOOP
	{
    r = sqrt( pow(p->XP[IP]-p->F58_1,2.0)+pow(p->YP[JP]-p->F58_2,2.0)+pow(p->ZP[KP]-p->F58_3,2.0));

	if(r<=p->F58_4)
	a->phi(i,j,k)=1.0;
	}
}

if(p->F59_r>0.0)
{
    
	LOOP
	{
    r = sqrt( pow(p->XP[IP]-p->F59_xm,2.0)+pow(p->YP[JP]-p->F59_ym,2.0));

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
        istart = p->posc_i(p->F70_xs[qn]);
        iend = p->posc_i(p->F70_xe[qn]);
        
        jstart = p->posc_j(p->F70_ys[qn]);
        jend = p->posc_j(p->F70_ye[qn]);
        
        kstart = p->posc_k(p->F70_zs[qn]);
        kend = p->posc_k(p->F70_ze[qn]);


        LOOP
        if(i>=istart && i<iend && j>=jstart && j<jend && k>=kstart && k<kend)
        a->phi(i,j,k)=1;
    }
	
	for(qn=0;qn<p->F71;++qn)
    {
        istart = p->posc_i(p->F71_xs[qn]);
        iend = p->posc_i(p->F71_xe[qn]);
        
        jstart = p->posc_j(p->F71_ys[qn]);
        jend = p->posc_j(p->F71_ye[qn]);
        
        kstart = p->posc_k(p->F71_zs[qn]);
        kend = p->posc_k(p->F71_ze[qn]);


        LOOP
        if(i>=istart && i<iend && j>=jstart && j<jend && k>=kstart && k<kend)
        a->phi(i,j,k)=-1;
    }
	
	for(qn=0;qn<p->F72;++qn)
    {
		istart = p->posc_i(p->F72_xs[qn]);
        iend = p->posc_i(p->F72_xe[qn]);
        
        jstart = p->posc_j(p->F72_ys[qn]);
        jend = p->posc_j(p->F72_ye[qn]);

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
	epsi = (1.6/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);
        
	dx = (a->phi(i+1,j,k)-a->phi(i-1,j,k))/(p->DXP[IM1]+p->DXP[IP]);
	dy = (a->phi(i,j+1,k)-a->phi(i,j-1,k))/(p->DYP[JM1]+p->DYP[JP]);
	dz = (a->phi(i,j,k+1)-a->phi(i,j,k-1))/(p->DZP[KM1]+p->DZP[KP]);
	
	dnorm = sqrt(p->DXN[IP]*p->DXN[IP] + p->DYN[JP]*p->DYN[JP] + p->DZN[KP]*p->DZN[KP]);
	
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
