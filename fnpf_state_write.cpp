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

#include"fnpf_state.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include<iostream>
#include<fstream>
#include<sys/stat.h>
#include<sys/types.h>

void fnpf_state::write(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
 
    // Open File
	int num=0;

    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = printcount;
	
	if(p->P40==1)
	num=0;
    
    //header
    header_ini(p,c,pgc);
    
    // result file
    filename(p,c,pgc,num);
	 
	ofstream result;
	result.open(name, ios::binary);
	
    // origin_xyz
    ddn=p->originx;
    result.write((char*)&ddn, sizeof (double));
    
    ddn=p->originy;
    result.write((char*)&ddn, sizeof (double));
    
    ddn=p->originz;
    result.write((char*)&ddn, sizeof (double));
    
    // origin_ijk
    iin=p->origin_i;
    result.write((char*)&iin, sizeof (int));
    
    iin=p->origin_j;
    result.write((char*)&iin, sizeof (int));
    
    iin=p->origin_k;
    result.write((char*)&iin, sizeof (int));
    
    // Nx,Ny,Nz
    iin=p->knox;
    result.write((char*)&iin, sizeof (int));
    
    iin=p->knoy;
    result.write((char*)&iin, sizeof (int));
    
    iin=p->knoz;
    result.write((char*)&iin, sizeof (int));
    
    

    
    ILOOP
    {
    ffn=p->XP[IP];
    result.write((char*)&ffn, sizeof (float));
    } 
    
    JLOOP
    {
    ffn=p->YP[JP];
    result.write((char*)&ffn, sizeof (float));
    } 
    
    KLOOP
    {
    ffn=p->ZP[KP];
    result.write((char*)&ffn, sizeof (float));
    } 
    
    SLICELOOP4
    {
    ffn=c->eta(i,j);
    result.write((char*)&ffn, sizeof (float));
    } 
    
    LOOP
    {
    ffn=c->u(i,j,k);
    result.write((char*)&ffn, sizeof (float));
    } 

	LOOP
    {
    ffn=c->v(i,j,k);
    result.write((char*)&ffn, sizeof (float));
    } 

	LOOP
    {
    ffn=c->w(i,j,k);
    result.write((char*)&ffn, sizeof (float));
    } 
	
	
	
	
	result.close();
	
	++printcount;
}