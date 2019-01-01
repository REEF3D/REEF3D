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

#include"pjm.h"
#include"lexer.h"
#include"fdm.h"

void pjm::debug(lexer* p,fdm* a)
{

    double x,y,z;

    char name[100];
    sprintf(name,"density_debug-%d.dat",p->mpirank+1);
	ofstream result;
	result.open(name);
    count=0;
	
	
	LOOP
	if(i==5)
	{
	result<<k<<" "<<roface(p,a,0,0,-1)<<endl;
		
		
	}

    

}
