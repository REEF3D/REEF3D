/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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

#include<iostream>
#include<fstream>
#include"6DOF_gc.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_gc::print_stl(lexer *p, fdm *a, ghostcell *pgc)
{
	int num=0;
	
	if(p->P15==1)
    num = p->printcount_sixdof;

    if(p->P15==2)
    num = p->count;
	
	if(num<0)
	num=0;
	
if(p->mpirank==0 && (((p->count%p->P20==0) && p->P30<0.0)  || (p->simtime>printtime && p->P30>0.0)   || p->count==0))
{
	printtime+=p->P30;
	
	if(p->P14==1)
	{
		if(num<10)
		sprintf(name,"./REEF3D_CFD_6DOF_STL/REEF3D-6DOF-00000%i.stl",num);

		if(num<100&&num>9)
		sprintf(name,"./REEF3D_CFD_6DOF_STL/REEF3D-6DOF-0000%i.stl",num);

		if(num<1000&&num>99)
		sprintf(name,"./REEF3D_CFD_6DOF_STL/REEF3D-6DOF-000%i.stl",num);

		if(num<10000&&num>999)
		sprintf(name,"./REEF3D_CFD_6DOF_STL/REEF3D-6DOF-00%i.stl",num);

		if(num<100000&&num>9999)
		sprintf(name,"./REEF3D_CFD_6DOF_STL/REEF3D-6DOF-0%i.stl",num);

		if(num>99999)
		sprintf(name,"./REEF3D_CFD_6DOF_STL/REEF3D-6DOF-%i.stl",num);
	}

	ofstream result;
	result.open(name, ios::binary);

	
	result<<"solid"<<" "<<"ascii"<<endl;
	
	double zero=0.0;
	
	for(n=0; n<tricount; ++n)
	{
	result<<" facet normal "<<zero<<" "<<zero<<" "<<zero<<endl;
	result<<"  outer loop"<<endl;
	result<<"   vertex "<<tri_x[n][0]<<" "<<tri_y[n][0]<<" "<<tri_z[n][0]<<endl;
	result<<"   vertex "<<tri_x[n][1]<<" "<<tri_y[n][1]<<" "<<tri_z[n][1]<<endl;
	result<<"   vertex "<<tri_x[n][2]<<" "<<tri_y[n][2]<<" "<<tri_z[n][2]<<endl;
	result<<"  endloop"<<endl;
	result<<" endfacet"<<endl;
	}
	
	result<<"endsolid"<<endl;
	

	result.close();

	++p->printcount_sixdof;	
}
}
