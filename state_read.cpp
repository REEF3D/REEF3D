/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"state.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"
#include<iostream>
#include<fstream>
#include<sys/stat.h>
#include<sys/types.h>

/*
  * u,v,w
  * press
  * phi
  * kin
  * eps
  * eddyv
  * topo
  * cbed
  * conc
 */

void state::read(lexer *p, fdm *a, ghostcell *pgc, turbulence *pturb)
{
    // Open File
	filename(p,a,pgc,p->I41);
	
	
	ifstream result;
	result.open(name, ios::binary);
	
    result.read((char*)&iin, sizeof (int));
	p->count=p->count_statestart=iin;
	
    result.read((char*)&iin, sizeof (int));
	p->printcount=iin;
	
    result.read((char*)&ddn, sizeof (double));
	p->simtime=ddn;
    
    result.read((char*)&ddn, sizeof (double));
	p->printtime=ddn;
    
    result.read((char*)&ddn, sizeof (double));
	p->sedprinttime=ddn;
    
    result.read((char*)&ddn, sizeof (double));
	p->fsfprinttime=ddn;
    
    result.read((char*)&ddn, sizeof (double));
	p->probeprinttime=ddn;
    
    result.read((char*)&ddn, sizeof (double));
	p->stateprinttime=ddn;
    
    ULOOP
    {
    result.read((char*)&ffn, sizeof (float));
    a->u(i,j,k)=double(ffn);
    }
	
	VLOOP
    {
    result.read((char*)&ffn, sizeof (float));
    a->v(i,j,k)=double(ffn);
    }
	
	WLOOP
    {
    result.read((char*)&ffn, sizeof (float));
    a->w(i,j,k)=double(ffn);
    }
	
	LOOP
    {
    result.read((char*)&ffn, sizeof (float));
    a->press(i,j,k)=double(ffn);
    }
	
	LOOP
    {
    result.read((char*)&ffn, sizeof (float));
    a->phi(i,j,k)=double(ffn);
    }
	
	LOOP
    {
    result.read((char*)&ffn, sizeof (float));
    pturb->kinget(i,j,k,ffn);
    }
	
	LOOP
    {
    result.read((char*)&ffn, sizeof (float));
    pturb->epsget(i,j,k,ffn);
    }
	
	LOOP
    {
    result.read((char*)&ffn, sizeof (float));
    a->eddyv(i,j,k)=double(ffn);
    }
	
	ALOOP
    {
    result.read((char*)&ffn, sizeof (float));
    a->topo(i,j,k)=double(ffn);
    }
	
	SLICELOOP4
    {
    result.read((char*)&ffn, sizeof (float));
    a->bedload(i,j)=double(ffn);
    }
	
	LOOP
    {
    result.read((char*)&ffn, sizeof (float));
    a->conc(i,j,k)=double(ffn);
    }
	
	int gcval_press, gcval_phi, gcval_topo;
	
	if(p->B76==0)
    gcval_press=40;  

    if(p->B76==1)
    gcval_press=41;

    if(p->B76==2)
    gcval_press=42;

    if(p->B76==3)
    gcval_press=43;
	
	if(p->B76==4) 
    gcval_press=44;
	
	if(p->B76==5) 
    gcval_press=45;
	
	
	if(p->F50==1)
	gcval_phi=51;

	if(p->F50==2)
	gcval_phi=52;

	if(p->F50==3)
	gcval_phi=53;

	if(p->F50==4)
	gcval_phi=54;
	
	
	if(p->S50==1)
	gcval_topo=151;

	if(p->S50==2)
	gcval_topo=152;

	if(p->S50==3)
	gcval_topo=153;
	
	if(p->S50==4)
	gcval_topo=154;
	
	
	pgc->start1(p,a->u,10);
	pgc->start2(p,a->v,11);
	pgc->start3(p,a->w,12);
    pgc->start4(p,a->press,gcval_press);
	pgc->start4(p,a->phi,gcval_phi);
	pturb->gcupdate(p,a,pgc);
	pgc->start4(p,a->eddyv,24);
	pgc->start4a(p,a->topo,gcval_topo);
	pgc->start4(p,a->conc,40);
	
	result.close();
}








