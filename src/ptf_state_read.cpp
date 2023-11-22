
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

#include"ptf_state.h"
#include"lexer.h"
#include"fdm_ptf.h"
#include"ghostcell.h"
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

void ptf_state::read(lexer *p, fdm_ptf *e, ghostcell *pgc)
{
    // Open File
	filename(p,e,pgc,p->I41);
	
	
	ifstream result;
	result.open(name, ios::binary);
	
    result.read((char*)&iin, sizeof (int));
	p->count=p->count_statestart=iin;
	
    result.read((char*)&iin, sizeof (int));
	p->printcount=iin-1;
    
    p->printcount = MAX(p->printcount,0);
	
    result.read((char*)&ddn, sizeof (double));
	p->simtime=ddn;
    
    result.read((char*)&ddn, sizeof (double));
	p->printtime=ddn;
    
    /*
    result.read((char*)&ddn, sizeof (double));
	p->sedprinttime=ddn;
    */
    result.read((char*)&ddn, sizeof (double));
	p->fsfprinttime=ddn;
    
    result.read((char*)&ddn, sizeof (double));
	p->probeprinttime=ddn;
    
    result.read((char*)&ddn, sizeof (double));
	p->stateprinttime=ddn;
    
    
    ALOOP
    {
    result.read((char*)&ffn, sizeof (float));
    e->topo(i,j,k)=double(ffn);
    }
    
    pgc->start4a(p,e->topo,150);
    
    // topoupdate
    //pgc->topo_update(p,e);
    
    ULOOP
    {
    result.read((char*)&ffn, sizeof (float));
    e->u(i,j,k)=double(ffn);
    }
	
	VLOOP
    {
    result.read((char*)&ffn, sizeof (float));
    e->v(i,j,k)=double(ffn);
    }
	
	WLOOP
    {
    result.read((char*)&ffn, sizeof (float));
    e->w(i,j,k)=double(ffn);
    }
	
	LOOP
    {
    result.read((char*)&ffn, sizeof (float));
    e->press(i,j,k)=double(ffn);
    }
	
	LOOP
    {
    result.read((char*)&ffn, sizeof (float));
    e->phi(i,j,k)=double(ffn);
    }
	/*
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
	*/
	LOOP
    {
    result.read((char*)&ffn, sizeof (float));
    e->eddyv(i,j,k)=double(ffn);
    }
	/*
	SLICELOOP4
    {
    result.read((char*)&ffn, sizeof (float));
    psed->qbeget(i,j,double(ffn));
    }
	*/
	LOOP
    {
    result.read((char*)&ffn, sizeof (float));
    e->conc(i,j,k)=double(ffn);
    }
	
	int gcval_press, gcval_phi, gcval_topo;
	
    gcval_press=40;  
	
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
	
	
	pgc->start1(p,e->u,10);
	pgc->start2(p,e->v,11);
	pgc->start3(p,e->w,12);
    pgc->start4(p,e->press,gcval_press);
	pgc->start4(p,e->phi,gcval_phi);
	//pturb->gcupdate(p,e,pgc);
	pgc->start4(p,e->eddyv,24);
	pgc->start4a(p,e->topo,gcval_topo);
	pgc->start4(p,e->conc,40);
	
	result.close();
}








