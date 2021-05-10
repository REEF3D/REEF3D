/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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
    // header file
    if(ini_token==0)
    {
    if(p->mpirank==0)
    mainheader_ini(p,c,pgc);
    
    header(p,c,pgc);
    
    ini_token=1;
    }
    
    if(p->mpirank==0)
    mainheader(p,c,pgc);
    
    
    
    
    // Open File
	int num=0;

    if(p->P15>=1)
    num = printcount;

    
    // result file
    filename(p,c,pgc,num);
	 
    if(p->mpirank==0){ 
	ofstream result;
	result.open(name, ios::binary);
    
     
    // head section
    iin=file_version;
    result.write((char*)&iin, sizeof (int));
    
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
    
    // result section
    SLICELOOP4
    {
    ffn=float(c->eta(i,j));
    result.write((char*)&ffn, sizeof (float));
    } 
    
    SLICELOOP4
    {
    ffn=float(c->Fifsf(i,j));
    result.write((char*)&ffn, sizeof (float));
    } 
    
    FLOOP
    {
    ffn=float(c->U[FIJK]);
    result.write((char*)&ffn, sizeof (float));
    } 

	FLOOP
    {
    ffn=float(c->V[FIJK]);
    result.write((char*)&ffn, sizeof (float));
    } 

	FLOOP
    {
    ffn=float(c->W[FIJK]);
    result.write((char*)&ffn, sizeof (float));
    } 
    
    if(p->P44==1)
    FLOOP
    {
    ffn=float(c->Fi[FIJK]);
    result.write((char*)&ffn, sizeof (float));
    } 
	
	
	result.close();}
	
	++printcount;
}
