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


void ptf_state::write_result(lexer *p, fdm_ptf *e, ghostcell *pgc)
{
    // Open File
	int num=0;

    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = printcount;
	
	if(p->P40==1)
	num=0;
    
    filename(p,e,pgc,num);
	 
	ofstream result;
	result.open(name, ios::binary);
	
	iin=p->count;
    result.write((char*)&iin, sizeof (int));
	
	iin=p->printcount;
    result.write((char*)&iin, sizeof (int));
	
	ddn=p->simtime;
    result.write((char*)&ddn, sizeof (double));
    
    ddn=p->printtime;
    result.write((char*)&ddn, sizeof (double));
    /*
    ddn=p->sedprinttime;
    result.write((char*)&ddn, sizeof (double));
    */
    ddn=p->fsfprinttime;
    result.write((char*)&ddn, sizeof (double));
    
    ddn=p->probeprinttime;
    result.write((char*)&ddn, sizeof (double));
    
    ddn=p->stateprinttime;
    result.write((char*)&ddn, sizeof (double));   
    
    ALOOP
    {
    ffn=e->topo(i,j,k);
    result.write((char*)&ffn, sizeof (float));
    } 
    
    ULOOP
    {
    ffn=e->u(i,j,k);
    result.write((char*)&ffn, sizeof (float));
    } 

	VLOOP
    {
    ffn=e->v(i,j,k);
    result.write((char*)&ffn, sizeof (float));
    } 
	
	WLOOP
    {
    ffn=e->w(i,j,k);
    result.write((char*)&ffn, sizeof (float));
    } 
	
	LOOP
    {
    ffn=e->press(i,j,k);
    result.write((char*)&ffn, sizeof (float));
    } 
	
	LOOP
    {
    ffn=e->phi(i,j,k);
    result.write((char*)&ffn, sizeof (float));
    } 
	/*
	LOOP
    {
    ffn=pturb->kinval(i,j,k);
    result.write((char*)&ffn, sizeof (float));
    } 
	
	LOOP
    {
    ffn=pturb->epsval(i,j,k);
    result.write((char*)&ffn, sizeof (float));
    } 
	*/
	LOOP
    {
    ffn=e->eddyv(i,j,k);
    result.write((char*)&ffn, sizeof (float));
    } 
    /*
	SLICELOOP4
    {
    ffn=psed->qbeval(i,j);
    result.write((char*)&ffn, sizeof (float));
    } 
	*/
	LOOP
    {
    ffn=e->conc(i,j,k);
    result.write((char*)&ffn, sizeof (float));
    } 
	
	
	
	result.close();
	
	++printcount;
}