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

fnpf_state::fnpf_state(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{	
	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_FNPF_State",0777);
	
	printcount=0;
    
    ini_token=0;
    
    is_flag=ie_flag=js_flag=je_flag=0;
    
    if(p->P43==1)
    {
        if(p->P43_xs>=p->originx && p->P43_xs<p->endx)
        {
        is = p->posc_i(p->P43_xs);
        is_flag=1;
        }
        
        
        if(p->P43_xe>=p->originx && p->P43_xe<p->endx)
        {
        ie = p->posc_i(p->P43_xe);
        ie_flag=1;
        }
        
        if(p->P43_ys>=p->originy && p->P43_ys<p->endy)
        {
        js = p->posc_j(p->P43_ys);
        js_flag=1;
        }
        
        
        if(p->P43_ye>=p->originy && p->P43_ye<p->endy)
        {
        je = p->posc_j(p->P43_ye);
        je_flag=1;
        }
        
        
    
        
    }
}

fnpf_state::~fnpf_state()
{
}
