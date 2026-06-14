/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"vrans_nhflow_f.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

void vrans_nhflow_f::print_force(lexer* p, fdm_nhf *d, ghostcell *pgc)
{
    // write to surf file
    fout<<p->simtime<<"\t";
    fout<<Fx<<" \t ";
    fout<<Fy<<" \t ";
	fout<<Fz;

    fout<<endl;
}

void vrans_nhflow_f::print_force_ini(lexer* p, fdm_nhf *d, ghostcell *pgc)
{
    // Create Folder
	if(p->mpirank==0)
	mkdir("./REEF3D_NHFLOW_VRANS_Force",0777);
	
    if(p->mpirank==0)
    {
    // open force surf file
	sprintf(name,"./REEF3D_NHFLOW_VRANS_Force/REEF3D_VRANS_NHFLOW_Force.dat");
	
	fout.open(name);

    
    fout<<"time \t Fx \t Fy \t Fz ";

    fout<<endl;
	}
}
