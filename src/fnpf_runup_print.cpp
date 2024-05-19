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

#include"fnpf_runup.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

void fnpf_runup::print_ini(lexer* p, fdm_fnpf *c, ghostcell *pgc)
{
    // Create Folder
	if(p->mpirank==0)
	mkdir("./REEF3D_FNPF_Runup",0777);
	
    if(p->mpirank==0)
    {
        // open fnpf_runup file
        sprintf(name,"./REEF3D_FNPF_Runup/REEF3D_Runup-%i.dat",ID+1);
        
        fout.open(name);


        fout<<p->P140_x[ID]<<" \t "<<p->P140_y[ID]<<" \t "<<endl;
        fout<<endl<<endl;
     
        fout<<"it \t time \t Fx \t Fy ";

        fout<<endl;
	}
}

void fnpf_runup::print_fnpf_runup(lexer* p, fdm_fnpf *c, ghostcell *pgc)
{
    // write to runup file
    fout<<p->count<<" \t "<<setprecision(9)<<p->simtime<<" \t "<<R1<<" \t "<<R2<<" \t "<<R3<<" \t "<<R4<<" \t "<<R5<<" \t "<<R6<<endl;
}
