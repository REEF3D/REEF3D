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

#include"directreini.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"solver.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"picard_f.h"
#include"picard_void.h"
#include"reini_RK3.h"
#include<sys/stat.h>
#include<sys/types.h>

directreini::directreini(lexer* p, fdm *a):gradient(p),vertice(p), nodeflag(p),d0(p),wallf(p),epsi(0.6*p->DXM),zero(0.0)
{
	if(p->F50==1)
	gcval_phi=51;

	if(p->F50==2)
	gcval_phi=52;

	if(p->F50==3)
	gcval_phi=53;

	if(p->F50==4)
	gcval_phi=54;

    gcval_iniphi=50;

    if(p->F61>1.0e-20 || p->F62>1.0e-20)
    gcval_iniphi=gcval_phi;

	gcval_ro=1;

	if(p->F46==1)
	ppicard = new picard_f(p);

	if(p->F46!=1)
	ppicard = new picard_void(p);

	ppreini = new reini_RK3(p,1);
	
	p->F49=0;
	p->F44=2;
	
	dT = p->F43*p->DXM;
	dV = pow(p->DXM,3.0);
	
}

directreini::~directreini()
{
}

void directreini::start(fdm* a,lexer* p,field& b, ghostcell* pgc,ioflow* pflow)
{
    starttime=pgc->timer();
	
	LOOP
	d0(i,j,k)=b(i,j,k);
	pgc->start4(p,d0,gcval_phi);

    ppicard->volcalc(p,a,pgc,a->phi);
    pgc->start4(p,b,gcval_phi);
	
	
//---------------
// Algorithm
    if(p->count>0)
    {
	triangulation(p,a,b,nodeflag,vertice);
	reconstruct(p,a,b,nodeflag,vertice);
	reini(p,a,pgc,b,nodeflag,vertice);
	constraint(p,a,pgc,b);
	//correction(p,a,pgc,b);
    }
	
	ppreini->start(a,p,b,pgc,pflow);
    //debug(p,a);

    if(p->count>0)
    finalize(p,a);
//---------------

	ppicard->correct_ls(p,a,pgc,a->phi);

	p->reinitime=pgc->timer()-starttime;
}

void directreini::startV(fdm* a,lexer* p,vec &f, ghostcell* pgc,ioflow* pflow)
{ 
    
}


