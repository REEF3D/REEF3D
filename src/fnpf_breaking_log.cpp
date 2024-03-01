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
#include"fnpf_breaking_log.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include<iostream>
#include<fstream>
#include<sys/stat.h>
#include<sys/types.h>

fnpf_breaking_log::fnpf_breaking_log(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{	
	// Create Folder
	mkdir("./REEF3D_FNPF_Breaking_Log",0777);
	
	// result file
    filename(p,c,pgc);
	
	result.open(name);
}

fnpf_breaking_log::~fnpf_breaking_log()
{
    result.close();
}

void fnpf_breaking_log::write(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    
    // result section
    SLICELOOP4
    if(c->breaklog(i,j)>0)
    {
    result<<p->simtime<<" "<<p->XP[IP]<<" "<<p->YP[JP]<<endl;
    //cout<<p->simtime<<" "<<p->XP[IP]<<" "<<p->YP[JP]<<endl;
    } 


}
