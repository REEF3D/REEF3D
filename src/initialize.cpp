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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"initialize.h"
#include"fdm.h"
#include"lexer.h"
#include"ghostcell.h"

initialize::initialize(lexer* p):smallnum(1.0e-20)
{
}

initialize::~initialize()
{
}

void initialize::start(fdm* a, lexer* p, ghostcell* pgc)
{
	deltax=p->DXM;
	
	inifdm(a,p,pgc);
	nodecalc(a,p);
	maxcoor(a,p,pgc);
	paraini(p,a,pgc);
    
    
    p->phimean=p->F56;
    
    if(p->F60>-1.0e20)
    p->phimean=p->F60;
        

	
	if(p->F40>0)
	iniphi(a,p,pgc);

	if(p->F80>0 && p->F80<4)
	inivof(a,p,pgc);
    
	if(p->F80==4)
	inivofPLIC(a,p,pgc);  

	if((p->F70>0 || p->F71>0 ||p->F72>0) && p->F40>0)
	iniphi_box(p,a,pgc);

	if(p->F70>0 && p->F80>0)
	inivof_box(p,a,pgc);
    	
	//iniphi_surfarea(p,a,pgc);

	if(p->S10>0)
	topoini(p,a,pgc);
	
	pgc->flagbase(p,a);
}

void initialize::inifdm(fdm* a, lexer* p, ghostcell* pgc)
{	
		ULOOP
		a->u(i,j,k)=0.0;

		VLOOP
		a->v(i,j,k)=0.0;

		WLOOP
        a->w(i,j,k)=0.0;

	LOOP
	{
		a->F(i,j,k)=0.0;
		a->G(i,j,k)=0.0;
		a->H(i,j,k)=0.0;

		a->press(i,j,k)=p->I55;
        
        a->Fi(i,j,k)=0.0;

		a->ro(i,j,k)=p->W1;
		a->visc(i,j,k)=p->W2;
		a->eddyv(i,j,k)=0.0;
		a->phi(i,j,k)=1.0;

		a->conc(i,j,k)=0.0;
	}

	ALOOP
    {
	a->fb(i,j,k)=1.0;
    a->topo(i,j,k)=1.0;
    a->porosity(i,j,k)=1.0;
    }

	pgc->start4(p,a->ro,1);
    pgc->start4(p,a->ro,1);
	pgc->start4(p,a->visc,1);
	pgc->start4(p,a->eddyv,1);
    pgc->start4a(p,a->porosity,1);
	pgc->start4a(p,a->press,1);
	pgc->start4a(p,a->fb,150);
	pgc->start4a(p,a->topo,150);
	pgc->gcparacox(p,a->fb,150);
    pgc->start4(p,a->phi,50);
}

void initialize::nodecalc(fdm* a, lexer* p)
{
	int count=0;
	p->pointnum=0;
	p->cellnum=0;
	i=0;
    
    // 3D
	TPLOOP
	{
	++count;
	++p->pointnum;
	a->nodeval(i,j,k)=count;
	}

	LOOP
	++p->cellnum;
    
    LOOP
    ++p->tpcellnum;
    
    // 2D
    count=0;
    TPSLICELOOP
	{
	++count;
	++p->pointnum2D;
	a->nodeval2D(i,j)=count;
    }
}

void initialize::maxcoor(fdm* a, lexer* p, ghostcell* pgc)
{
p->maxlength=-1.0e9;
p->xcoormax=-1.0e9;
p->xcoormin=1.0e9;
p->ycoormax=-1.0e9;
p->ycoormin=1.0e9;
p->zcoormax=-1.0e9;
p->zcoormin=1.0e9;

    LOOP
    {
        p->xcoormax = MAX(p->xcoormax,p->XN[IP1]);
        p->xcoormin = MIN(p->xcoormin,p->XN[IP]);
        p->ycoormax = MAX(p->ycoormax,p->YN[JP1]);
        p->ycoormin = MIN(p->ycoormin,p->YN[JP]);
        p->zcoormax = MAX(p->zcoormax,p->ZN[KP1]);
        p->zcoormin = MIN(p->zcoormin,p->ZN[KP]);
     }

     p->maxlength=MAX(p->maxlength,p->xcoormax-p->xcoormin);
     p->maxlength=MAX(p->maxlength,p->ycoormax-p->ycoormin);
     p->maxlength=MAX(p->maxlength,p->zcoormax-p->zcoormin);

     p->maxlength=pgc->globalmax(p->maxlength);
	 
	 p->xcoormax=pgc->globalmax(p->xcoormax);
	 p->ycoormax=pgc->globalmax(p->ycoormax);
	 p->zcoormax=pgc->globalmax(p->zcoormax);
	 
	 p->xcoormin=pgc->globalmin(p->xcoormin);
	 p->ycoormin=pgc->globalmin(p->ycoormin);
	 p->zcoormin=pgc->globalmin(p->zcoormin);
	 
	 if(p->F42>=0.0)
	 p->maxlength = p->F42;
}

int initialize::conv(double a)
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

