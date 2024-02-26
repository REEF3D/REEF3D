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
#include"fnpf_ini.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"

void fnpf_ini::fnpf_restart(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{

// Open File
    // read mainheader, find file_type
    fnpf_restart_mainheader(p,c,pgc);
    
    // if single file
    if(file_type==1)
	filename(p,c,pgc,p->I41);
    
    // if contiuous file
    if(file_type==2)
    {
        sprintf(name,"./REEF3D_FNPF_STATE/REEF3D_FNPF-State-%06i.r3d",p->mpirank+1);
    }
	
    // open file
	result.open(name, ios::binary);
	
    // read single state file
    if(file_type==1) 
    fnpf_restart_read(p,c,pgc);
    
    // read continuous state file
    if(file_type==2)
    do
    {
        fnpf_restart_read(p,c,pgc);
         
    }while(p->count_statestart<p->I41);

	result.close();
    
    
    // ----------------------
    // finish: ghostell update
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
	

}

