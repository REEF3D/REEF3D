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

#include"initialise_ptf.h"
#include"fdm_ptf.h"
#include"lexer.h"
#include"ghostcell.h"

initialise_ptf::initialise_ptf(lexer* p):smallnum(1.0e-20)
{
}

initialise_ptf::~initialise_ptf()
{
}

void initialise_ptf::start(fdm_ptf *e, lexer* p, ghostcell* pgc)
{
	deltax=p->DXM;
	
	inifdm_ptf(e,p,pgc);
	nodecalc(e,p);
	maxcoor(e,p,pgc);
	paraini(p,e,pgc);
    
    
    p->phimean=p->F56;
    
    if(p->F60>-1.0e20)
    p->phimean=p->F60;
        

	
	if(p->F40>0)
	iniphi(e,p,pgc);

	if((p->F70>0 || p->F71>0 ||p->F72>0) && p->F40>0)
	iniphi_box(p,e,pgc);
    	
	//iniphi_surfarea(p,e,pgc);

	if(p->S10>0)
	topoini(p,e,pgc);
	
	pgc->flagbase_ptf(p,e);
}

void initialise_ptf::inifdm_ptf(fdm_ptf *e, lexer* p, ghostcell* pgc)
{	
		ULOOP
		e->u(i,j,k)=0.0;

		VLOOP
		e->v(i,j,k)=0.0;

		WLOOP
        e->w(i,j,k)=0.0;

	LOOP
	{
		e->F(i,j,k)=0.0;
		e->G(i,j,k)=0.0;
		e->H(i,j,k)=0.0;

		e->press(i,j,k)=p->I55;
        
        e->Fi(i,j,k)=0.0;

		e->ro(i,j,k)=p->W1;
		e->visc(i,j,k)=p->W2;
		e->eddyv(i,j,k)=0.0;
		e->phi(i,j,k)=1.0;

		e->conc(i,j,k)=0.0;
	}

	ALOOP
    {
	e->fb(i,j,k)=1.0;
    e->topo(i,j,k)=1.0;
    e->porosity(i,j,k)=1.0;
    }

	pgc->start4(p,e->ro,1);
    pgc->start4(p,e->ro,1);
	pgc->start4(p,e->visc,1);
	pgc->start4(p,e->eddyv,1);
    pgc->start4a(p,e->porosity,1);
	pgc->start4a(p,e->press,1);
	pgc->start4a(p,e->fb,150);
	pgc->start4a(p,e->topo,150);
	pgc->gcparacox(p,e->fb,150);
    pgc->start4(p,e->phi,50);
}

void initialise_ptf::nodecalc(fdm_ptf *e, lexer* p)
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
	e->nodeval(i,j,k)=count;
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
	e->nodeval2D(i,j)=count;
    }
}

void initialise_ptf::maxcoor(fdm_ptf *e, lexer* p, ghostcell* pgc)
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

int initialise_ptf::conv(double a)
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

