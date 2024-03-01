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
#include"gradient.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include<sys/stat.h>
#include<sys/types.h>

force_ale::force_ale(lexer* p, fdm_fnpf *c, ghostcell *pgc, int qn) : ID(qn){}

force_ale::~force_ale(){}

void force_ale::ini(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    force_aleprintcount=0;

    // Read cylinder force input - xc,yc,rc,cd,cm
    xc = p->P85_x[ID];
	yc = p->P85_y[ID];
    rc = p->P85_r[ID];
	cd = p->P85_cd[ID];
	cm = p->P85_cm[ID];

    // Open files
    print_ini(p,c,pgc);

    // Ini arrays
	p->Darray(un, p->knoz);
	//p->Darray(u2n, p->knoz);
	p->Darray(vn, p->knoz);

    // Ini eta
	etan=p->wd;
	//eta2n=p->wd;

    // Ini time
    //dtn=0;

    // Ini processor boundaries
	xstart = p->originx;
	ystart = p->originy;
	xend = p->endx;
	yend = p->endy;

}

void force_ale::start(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if (xc >= xstart && xc < xend && yc >= ystart && yc < yend) // cylinder in processor
    {
        i = p->posc_i(xc);
        j = p->posc_j(yc);

        // Calculate force
        force_ale_force(p,c,pgc);
    }
    else
    {
        Fx = Fy = 0.0;
    }

    // Sum up to distribute forces
    Fx = pgc->globalsum(Fx);
    Fy = pgc->globalsum(Fy);

    // Print
    if(p->mpirank==0)
    {
        print_force_ale(p,c,pgc);
    }
}
