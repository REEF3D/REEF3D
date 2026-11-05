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

#include<sys/stat.h>
#include<iostream>
#include<fstream>
#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_obj::print_parameter(lexer *p, ghostcell *pgc)
{
    print_parameter(p,pgc,p->simtime);
}

void sixdof_obj::print_parameter(lexer *p, ghostcell *pgc, double t)
{
    // t: time of the printed state. NHFLOW prints the body at the end of the step, before the
    // driver advances simtime, so it passes simtime + dt (the label of the fluid output)
	if(p->mpirank==0 && p->count%p->X19==0)
    {
        ofstream print;
        char str[1000];
        
        
        printpos<<t<<" \t "<<p->xg<<" \t "<<p->yg<<" \t "<<p->zg<<" \t "<<phi*(180/PI)<<" \t "<<theta*(180/PI)<<" \t "<<psi*(180/PI)<<endl;


        printvel<<t<<" \t "<<p->ufbi<<" \t "<<p->vfbi<<" \t "<<p->wfbi<<" \t "<<p->pfbi<<" \t "<<p->qfbi<<" \t "<<p->rfbi<<endl;

    }
    
    // load models (ship module)
    for(size_t ql=0; ql<pload.size(); ++ql)
    pload[ql]->print(p);
}
