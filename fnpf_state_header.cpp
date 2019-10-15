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

void fnpf_state::header_ini(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{

  // Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_FNPF_STATE",0777);
	
	hdout.open(name, ios::binary);
    
    // open file
	if(p->P14==0)
    hdout.open("REEF3D-FNPF_state_header.r3d");
	
	if(p->P14==1)
	hdout.open("./REEF3D_FNPF_STATE/REEF3D-FNPF_state_header.r3d");
    
    
    // ini write
    iin=p->M10;
    hdout.write((char*)&iin, sizeof (int));
    
    iin=p->gknox;
    hdout.write((char*)&iin, sizeof (int));
    
    iin=p->gknoy;
    hdout.write((char*)&iin, sizeof (int));
    
    iin=p->gknoz+1;
    hdout.write((char*)&iin, sizeof (int));
}


void fnpf_state::header(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    iin=p->count;
    hdout.write((char*)&iin, sizeof (int));
		
	ddn=p->simtime;
    hdout.write((char*)&ddn, sizeof (double));
}