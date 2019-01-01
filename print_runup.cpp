/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"print_runup.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

print_runup::print_runup(lexer *p, fdm* a, ghostcell *pgc)
{	
	
	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_Runup-FSF",0777);
	
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_Runup-MAX",0777);
	
	cone_ini(p,a,pgc);
}

print_runup::~print_runup()
{
}


void print_runup::start(lexer *p, fdm* a, ghostcell *pgc)
{
	cone(p,a,pgc);
	print(p,a,pgc);
}














