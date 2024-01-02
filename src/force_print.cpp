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

#include"force.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

void force::print_force(lexer* p, fdm *a, ghostcell *pgc)
{
    // write to surf file

    fout<<p->count<<"\t";
    fout<<setprecision(9)<<p->simtime<<"\t";
    fout<<Fx<<" \t ";
    fout<<Fy<<" \t ";
	fout<<Fz;
    

    fout<<endl;
}

void force::print_ini(lexer* p, fdm *a, ghostcell *pgc)
{
    // Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_CFD_Force",0777);
	
    if(p->mpirank==0)
    {
    // open force surf file
	if(p->P14==0)
	sprintf(name,"REEF3D_CFD_Force-%i.dat",ID+1);
    
	if(p->P14==1)
	sprintf(name,"./REEF3D_CFD_Force/REEF3D_CFD_Force-%i.dat",ID+1);
	
	fout.open(name);

    fout<<"x_start xend     y_start y_end     z_start z_end"<<endl;

    fout<<p->P81_xs[ID]<<" "<<p->P81_xe[ID]<<" . "<<p->P81_ys[ID]<<" "<<p->P81_ye[ID]<<" . "<<p->P81_zs[ID]<<" "<<p->P81_ze[ID]<<endl;
    fout<<endl<<endl;
    
 
    fout<<"it \t time \t Fx \t Fy \t Fz ";
    

    fout<<endl;
	}

    
}
