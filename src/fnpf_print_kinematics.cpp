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

#include"fnpf_print_kinematics.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include<sys/stat.h>
#include<sys/types.h>

fnpf_print_kinematics::fnpf_print_kinematics(lexer* p, fdm_fnpf *c, ghostcell *pgc, int qn) : ID(qn)
{
    // Create Folder
	if(p->mpirank==0)
	mkdir("./REEF3D_FNPF_Kinematics",0777);
}

fnpf_print_kinematics::~fnpf_print_kinematics()
{
}

void fnpf_print_kinematics::ini(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    fnpf_print_kinematicsprintcount=0;
    
    // Ini processor boundaries
	xstart = p->originx;
	ystart = p->originy;
	xend = p->endx;
	yend = p->endy;

    // Read cylinder force input - xc,yc,rc,cd,cm
    xc = p->P88_x[ID];
	yc = p->P88_y[ID];

    // Open files
    if (xc >= xstart && xc < xend && yc >= ystart && yc < yend) 
    print_ini(p,c,pgc);

    // Ini arrays
	p->Darray(un, p->knoz+1);
	p->Darray(vn, p->knoz+1);
    p->Darray(ax, p->knoz+1);
	p->Darray(ay, p->knoz+1);

    // Ini eta
	etan=p->wd;
}

void fnpf_print_kinematics::start(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if (xc >= xstart && xc < xend && yc >= ystart && yc < yend) 
    {
        i = p->posc_i(xc);
        j = p->posc_j(yc);

        // calculate kinematics
        kinematics_calc(p,c,pgc);
        
        // print kinematics
        print_kinematics(p,c,pgc);
    }
}
