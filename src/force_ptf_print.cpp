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
Author: Fabian Knoblauch
--------------------------------------------------------------------*/
#include"force_ptf.h"
#include"lexer.h"
#include"fdm_ptf.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

void force_ptf::print_force_ptf(lexer* p, fdm_ptf *e, ghostcell *pgc)
{
    // write to force file
    if(p->P87==0)
        fout<<p->count<<" \t "<<setprecision(9)<<p->simtime<<" \t "<<F_x_tot<<" \t "<<F_y_tot<<" \t "<<F_z_tot<<endl;
    if(p->P87==1)
        fout<<setprecision(9)<<p->simtime<<","<<F_x_tot<<","<<F_y_tot<<","<<F_z_tot<<endl;
}

void force_ptf::print_ini_ptf(lexer* p, fdm_ptf *e, ghostcell *pgc)
{
    // Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_PTF_Force_PTF",0777);
	
    if(p->mpirank==0)
    {
        // open force_ptf file
        if(p->P14==0)
        {
            if(p->P87==0)
                sprintf(name,"REEF3D_PTF_Force-%i.dat",ID+1);
            if(p->P87==1)
                sprintf(name,"REEF3D_PTF_Force-%i.csv",ID+1);
        }
        
        if(p->P14==1)
        {
            if(p->P87==0)
                sprintf(name,"./REEF3D_PTF_Force_PTF/REEF3D_PTF_Force.dat",ID+1);
            if(p->P87==1)
                sprintf(name,"./REEF3D_PTF_Force_PTF/REEF3D_PTF_Force.csv",ID+1);
        }
        
        if(p->P87==0)
        {
        fout.open(name);

        fout<<"it \t time \t Fx \t Fy \t Fz";

        fout<<endl;
        }
        
        if(p->P87==1)
        {
            fout.open(name);
            fout<<"time,Fx,Fy,Fz";
            fout<<endl;
        }
	}
}
