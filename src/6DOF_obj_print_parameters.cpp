/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
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
	if(p->mpirank==0 && p->count%p->X19==0)
    {
        ofstream print;
        char str[1000];
        
        
        printpos<<p->simtime<<" \t "<<p->xg<<" \t "<<p->yg<<" \t "<<p->zg<<" \t "<<phi*(180/PI)<<" \t "<<theta*(180/PI)<<" \t "<<psi*(180/PI)<<endl;


        printvel<<p->simtime<<" \t "<<p->ufbi<<" \t "<<p->vfbi<<" \t "<<p->wfbi<<" \t "<<p->pfbi<<" \t "<<p->qfbi<<" \t "<<p->rfbi<<endl;

    }
}
