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
--------------------------------------------------------------------*/

#include<sys/stat.h>
#include<iostream>
#include<fstream>
#include"6DOF_df_object.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_df_object::print_stl(lexer *p, fdm *a, ghostcell *pgc)
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
        
        char path[300];
        
        if(p->P14==1)
        {
            if(num<10)
            sprintf(path,"./REEF3D_CFD_6DOF_STL/REEF3D-6DOF-%i-00000%i.stl",n6DOF,num);

            if(num<100&&num>9)
            sprintf(path,"./REEF3D_CFD_6DOF_STL/REEF3D-6DOF-%i-0000%i.stl",n6DOF,num);

            if(num<1000&&num>99)
            sprintf(path,"./REEF3D_CFD_6DOF_STL/REEF3D-6DOF-%i-000%i.stl",n6DOF,num);

            if(num<10000&&num>999)
            sprintf(path,"./REEF3D_CFD_6DOF_STL/REEF3D-6DOF-%i-00%i.stl",n6DOF,num);

            if(num<100000&&num>9999)
            sprintf(path,"./REEF3D_CFD_6DOF_STL/REEF3D-6DOF-%i-0%i.stl",n6DOF,num);

            if(num>99999)
            sprintf(path,"./REEF3D_CFD_6DOF_STL/REEF3D-6DOF-%i-%i.stl",n6DOF,num);
        }

        ofstream result;
        result.open(path, ios::binary);

        
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


void sixdof_df_object::print_parameter(lexer *p, fdm *a, ghostcell *pgc)
{
	if(p->mpirank == 0 && p->count%p->X19==0)
    {
        ofstream print;
        char str[1000];
        
        if(p->P14==0)
        sprintf(str,"REEF3D_6DOF_position_%i.dat",n6DOF);
        if(p->P14==1)
        sprintf(str,"./REEF3D_CFD_6DOF/REEF3D_6DOF_position_%i.dat",n6DOF);
        
        print.open(str, std::ofstream::out | std::ofstream::app);
        print<<p->simtime<<" \t "<<p->xg<<" \t "<<p->yg<<" \t "<<p->zg<<" \t "<<phi*(180/PI)<<" \t "<<theta*(180/PI)<<" \t "<<psi*(180/PI)<<endl;
        print.close();
        
        
        if(p->P14==0)
        sprintf(str,"REEF3D_6DOF_velocity_%i.dat",n6DOF);
        if(p->P14==1)
        sprintf(str,"./REEF3D_CFD_6DOF/REEF3D_6DOF_velocity_%i.dat",n6DOF);
        
        print.open(str, std::ofstream::out | std::ofstream::app);
        print<<p->simtime<<" \t "<<p->ufbi<<" \t "<<p->vfbi<<" \t "<<p->wfbi<<" \t "<<p->pfbi<<" \t "<<p->qfbi<<" \t "<<p->rfbi<<endl;
        print.close();
    }
}
