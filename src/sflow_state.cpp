/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"sflow_state.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include<iostream>
#include<fstream>
#include<sys/stat.h>
#include<sys/types.h>

sflow_state::sflow_state(lexer *p, fdm2D *b, ghostcell *pgc, int state_restart)
{	
    restart=state_restart;
    
	// Create Folder
	if(p->mpirank==0 && restart==0)
	mkdir("./REEF3D_SFLOW_STATE",0777);
    
    if(p->mpirank==0 && restart==1)
	mkdir("./REEF3D_SFLOW_STATE_RESTART",0777);
	
	printcount=0;
    
    file_version=2;
    
    file_type=p->P45;
    
    ini_token=0;
    
    boundary(p,b,pgc,restart);
}

sflow_state::~sflow_state()
{
    result.close();
}

void sflow_state::write(lexer *p, fdm2D *c, ghostcell *pgc)
{
    // header file
    if(ini_token==0)
    {
    if(p->mpirank==0)
    ini_mainheader(p,c,pgc);
    
    if(flag==1) 
    write_header(p,c,pgc);
    
    ini_token=1;
    }
    
    if(p->mpirank==0)
    write_mainheader(p,c,pgc);
    
    
    // result file
    if(flag==1)
    write_result(p,c,pgc);
}

void sflow_state::write_single(lexer *p, fdm2D *c, ghostcell *pgc)
{
    
}

void sflow_state::write_contiuous(lexer *p, fdm2D *c, ghostcell *pgc)
{
    
}


void sflow_state::write_restart(lexer *p, fdm2D *c, ghostcell *pgc)
{
    
}
