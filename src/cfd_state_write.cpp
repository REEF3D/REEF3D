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

#include"cfd_state.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"
#include"sediment.h"
#include<iostream>
#include<fstream>
#include<sys/stat.h>
#include<sys/types.h>

void cfd_state::write_result(lexer *p, fdm *a, ghostcell *pgc, turbulence *pturb, sediment *psed)
{
    // Open File
	int num=0;

    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = printcount;
	
	if(p->P40==1)
	num=0;
    
    filename(p,a,pgc,num);
	 
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
    
    ddn=p->sedprinttime;
    result.write((char*)&ddn, sizeof (double));
    
    ddn=p->fsfprinttime;
    result.write((char*)&ddn, sizeof (double));
    
    ddn=p->probeprinttime;
    result.write((char*)&ddn, sizeof (double));
    
    ddn=p->stateprinttime;
    result.write((char*)&ddn, sizeof (double));   
    
    ALOOP
    {
    ffn=a->topo(i,j,k);
    result.write((char*)&ffn, sizeof (float));
    } 
    
    ULOOP
    {
    ffn=a->u(i,j,k);
    result.write((char*)&ffn, sizeof (float));
    } 

	VLOOP
    {
    ffn=a->v(i,j,k);
    result.write((char*)&ffn, sizeof (float));
    } 
	
	WLOOP
    {
    ffn=a->w(i,j,k);
    result.write((char*)&ffn, sizeof (float));
    } 
	
	LOOP
    {
    ffn=a->press(i,j,k);
    result.write((char*)&ffn, sizeof (float));
    } 
	
	LOOP
    {
    ffn=a->phi(i,j,k);
    result.write((char*)&ffn, sizeof (float));
    } 
	
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
	
	LOOP
    {
    ffn=a->eddyv(i,j,k);
    result.write((char*)&ffn, sizeof (float));
    } 

	SLICELOOP4
    {
    ffn=psed->qbeval(i,j);
    result.write((char*)&ffn, sizeof (float));
    } 
	
	LOOP
    {
    ffn=a->conc(i,j,k);
    result.write((char*)&ffn, sizeof (float));
    } 
	
	
	
	result.close();
	
	++printcount;
}