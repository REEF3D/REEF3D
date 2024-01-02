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

#include"fnpf_state.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include<iostream>
#include<fstream>
#include<sys/stat.h>
#include<sys/types.h>

void fnpf_state::write_result(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    // Open result file single file
    if(p->P45==1)
    {
    int num=0;

    if(p->P15>=1)
    num = printcount;
    
    filename_single(p,c,pgc,num);
	 
	result.open(name, ios::binary);
    }
     
    // read head section
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
    
    // read result section
    for(i=is;i<ie;++i)
    for(j=js;j<je;++j)
    PSLICECHECK4
    {
    ffn=float(c->eta(i,j));
    result.write((char*)&ffn, sizeof (float));
    } 
    
    for(i=is;i<ie;++i)
    for(j=js;j<je;++j)
    PSLICECHECK4
    {
    ffn=float(c->Fifsf(i,j));
    result.write((char*)&ffn, sizeof (float));
    } 
    
    for(i=is;i<ie;++i)
    for(j=js;j<je;++j)
    for(k=0; k<p->knoz+1; ++k)
    FPCHECK  
    {
    ffn=float(c->U[FIJK]);
    result.write((char*)&ffn, sizeof (float));
    } 

	for(i=is;i<ie;++i)
    for(j=js;j<je;++j)
    for(k=0; k<p->knoz+1; ++k)
    FPCHECK 
    {
    ffn=float(c->V[FIJK]);
    result.write((char*)&ffn, sizeof (float));
    } 

	for(i=is;i<ie;++i)
    for(j=js;j<je;++j)
    for(k=0; k<p->knoz+1; ++k)
    FPCHECK 
    {
    ffn=float(c->W[FIJK]);
    result.write((char*)&ffn, sizeof (float));
    } 
    
    if(p->P44==1)
    for(i=is;i<ie;++i)
    for(j=js;j<je;++j)
    for(k=0; k<p->knoz+1; ++k)
    FPCHECK 
    {
    ffn=float(c->Fi[FIJK]);
    result.write((char*)&ffn, sizeof (float));
    } 
	
	if(p->P45==1)
	result.close();
	
	++printcount;
}
