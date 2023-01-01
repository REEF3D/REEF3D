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

#include "export.h"
#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"
#include<iostream>
#include<fstream>
#include<sys/stat.h>
#include<sys/types.h>

exportfile::exportfile(lexer *p, fdm *a, ghostcell *pgc)
{	
	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_Export",0777);
	
	printcount=0;
    
    p->Darray(eta,p->knox,p->knoy);
}

exportfile::~exportfile()
{
}

void exportfile::start(lexer *p, fdm *a, ghostcell *pgc)
{
 
    preproc(p,a,pgc);
    
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
	
    // iteration count
	iin=p->count;
    result.write((char*)&iin, sizeof (int));
	
    // t
	ddn=p->simtime;
    result.write((char*)&ddn, sizeof (double));
    
    // dx
    ddn=p->DXM;
    result.write((char*)&ddn, sizeof (double));
    
    // xs,xe,ys,ye,zs,ze
    ddn=xs;
    result.write((char*)&ddn, sizeof (double));
    ddn=xe;
    result.write((char*)&ddn, sizeof (double));
    ddn=ys;
    result.write((char*)&ddn, sizeof (double));
    ddn=ye;
    result.write((char*)&ddn, sizeof (double));
    ddn=zs;
    result.write((char*)&ddn, sizeof (double));
    ddn=ze;
    result.write((char*)&ddn, sizeof (double));
    
    // Ni,Nj,Nk
    iin=p->knox;
    result.write((char*)&iin, sizeof (int));
    iin=p->knoy;
    result.write((char*)&iin, sizeof (int));
    iin=p->knoz;
    result.write((char*)&iin, sizeof (int));
 
    
    // u-velocity
    LOOP
    {
    ffn=0.5*(a->u(i,j,k)+a->u(i-1,j,k));
    result.write((char*)&ffn, sizeof (float));
    } 
    
    // v-velocity
	LOOP
    {
    ffn=0.5*(a->v(i,j,k)+a->v(i,j-1,k));
    result.write((char*)&ffn, sizeof (float));
    } 
	
    // w-velocity
	LOOP
    {
    ffn=0.5*(a->w(i,j,k)+a->w(i,j,k-1));
    result.write((char*)&ffn, sizeof (float));
    } 
	
    // phi
	LOOP
    {
    ffn=a->phi(i,j,k);
    result.write((char*)&ffn, sizeof (float));
    } 
    
    // eta : 2D 
	/*ILOOP
    JLOOP
    {
    ffn=eta[i][j];
    result.write((char*)&ffn, sizeof (float));
    } 
*/

	result.close();
	
	++printcount;
}
