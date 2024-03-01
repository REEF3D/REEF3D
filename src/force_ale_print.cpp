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

#include"force_ale.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

void force_ale::print_force_ale(lexer* p, fdm_fnpf *c, ghostcell *pgc)
{
    //cout<<"Fx"<<ID + 1<<": "<<Fx<<" Fy"<<ID + 1<<": "<<Fy<<endl;
    
    // write to force file
    fout<<p->count<<" \t "<<setprecision(9)<<p->simtime<<" \t "<<Fx<<" \t "<<Fy<<endl;
}

void force_ale::print_ini(lexer* p, fdm_fnpf *c, ghostcell *pgc)
{
    // Create Folder
	if(p->mpirank==0)
	mkdir("./REEF3D_FNPF_Force_ALE",0777);
	
    if(p->mpirank==0)
    {
        // open force_ale file
        sprintf(name,"./REEF3D_FNPF_Force_ALE/REEF3D_ALE_Force-%i.dat",ID+1);
        
        fout.open(name);

        fout<<"x \t y \t Cd \t Cm"<<endl;

        fout<<p->P85_x[ID]<<" \t "<<p->P85_y[ID]<<" \t "<<p->P85_cd[ID]<<" \t "<<p->P85_cm[ID] <<endl;
        fout<<endl<<endl;
     
        fout<<"it \t time \t Fx \t Fy ";

        fout<<endl;
	}
}
