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
#include"fnpf_ini.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"reini.h"
#include"convection.h"
#include"ioflow.h"

void fnpf_ini::restart(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{

    float ffn;
	int iin;
	double ddn;
	int printcount;

// Open File
	filename(p,c,pgc,p->I41);
	
	
	ifstream result;
	result.open(name, ios::binary);
	
    
    // head section
    result.read((char*)&iin, sizeof (int));
    
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
    
    // result section
    SLICELOOP4
    {
    result.read((char*)&ffn, sizeof (float));
    c->eta(i,j)=double(ffn);
    } 
    
    SLICELOOP4
    {
    result.read((char*)&ffn, sizeof (float));
    c->Fifsf(i,j)=double(ffn);
    } 
    
    FLOOP
    result.read((char*)&ffn, sizeof (float));
    
    FLOOP
    result.read((char*)&ffn, sizeof (float));
    
    FLOOP
    result.read((char*)&ffn, sizeof (float));
    
    if(p->I44==1)
    FLOOP
    result.read((char*)&ffn, sizeof (float));

	
    
    // ghostell update
	int gcval,gcval_u,gcval_v,gcval_w;
    int gcval_eta,gcval_fifsf;
    
    gcval=250;
    if(p->j_dir==0)
    gcval=150;
   
    gcval_u=10;
	gcval_v=11;
	gcval_w=12;
    
    // 3D
    gcval_eta = 55;
    gcval_fifsf = 60;
    
    // 2D
    if(p->j_dir==0)
    {
    gcval_eta = 155;
    gcval_fifsf = 160;
    }
	
    
    pgc->gcsl_start4(p,c->eta,gcval_eta);
    pgc->gcsl_start4(p,c->Fifsf,gcval_fifsf);
	

	result.close();
    
}

