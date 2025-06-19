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

#include"fnpf_runup.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"ioflow.h"

fnpf_runup::fnpf_runup(lexer* p, fdm_fnpf *c, ghostcell *pgc, int qn) : ID(qn)
{
}

fnpf_runup::~fnpf_runup()
{
}

void fnpf_runup::ini(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    fnpf_runupprintcount=0;

    // Read cylinder force input - xc,yc,rc
    xc = p->P140_x[ID];
	yc = p->P140_y[ID];
    rc = p->P141; //Edgar: Assign a varibale the radius of the cylinder

    // Open files
    print_ini(p,c,pgc);

    // Ini arrays
	//p->Darray(un, p->knoz);
	//p->Darray(vn, p->knoz);

    // Ini eta
	etan = p->wd;

    // Ini processor boundaries
	xstart = p->originx;
	ystart = p->originy;
	xend = p->endx;
	yend = p->endy;
}

void fnpf_runup::start(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    R1=R2=R3=R4=R5=R6=0.0;
    
    if (xc >= xstart && xc < xend && yc >= ystart && yc < yend) // cylinder in processor
    {
        i = p->posc_i(xc);
        j = p->posc_j(yc);
        
        cout<<"Run-up "<<i<<" "<<j<<endl;

        // Calculate runup
        fnpf_runup_calc(p,c,pgc);
    }

    // Sum up to distribute forces
    R1 = pgc->globalmax(R1);
    R2 = pgc->globalmax(R2);
    R3 = pgc->globalmax(R3);
    R4 = pgc->globalmax(R4);
    R5 = pgc->globalmax(R5);
    R6 = pgc->globalmax(R6);


    // Print
    if(p->mpirank==0)
    print_fnpf_runup(p,c,pgc);
}
