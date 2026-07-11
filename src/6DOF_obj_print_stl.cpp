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

void sixdof_obj::print_stl(lexer *p, ghostcell *pgc)
{
	int num=0;
    int printflag=0;
	
	if(p->P15==1)
    num = p->printcount_sixdof;

    if(p->P15==2)
    num = p->count;
	
	if(num<0)
	num=0;

    // CFD
    if(p->A10==2 || p->A10==5)
    {
        if(((p->count%p->P181==0) && p->P182<0.0) 
            || (p->simtime>printtime && p->P182>0.0) 
            || (p->count==0 && p->P185==0))
            printflag=true;

        if(p->P185>0)
        for(int qn=0; qn<p->P185; ++qn)
        if(p->simtime>printtime_wT[qn] 
            && p->simtime>=p->P185_ts[qn] 
            && p->simtime<=(p->P185_te[qn]+0.5*p->P185_dt[qn]))
        {
            printflag=true;

            printtime_wT[qn]+=p->P185_dt[qn];
        }
    }

    // CFD
    if(p->A10==6)
    {
        if(((p->count%p->P20==0) && p->P30<0.0) 
            || (p->simtime>printtime && p->P30>0.0) 
            || (p->count==0 && p->P35==0))
            printflag=true;

        if(p->P35>0)
        for(int qn=0; qn<p->P35; ++qn)
        if(p->simtime>printtime_wT[qn] 
            && p->simtime>=p->P35_ts[qn] 
            && p->simtime<=(p->P35_te[qn]+0.5*p->P35_dt[qn]))
        {
            printflag=true;

            printtime_wT[qn]+=p->P35_dt[qn];
        }
    }
    
    if(p->mpirank==0 && printflag==1)
    {
        if(p->A10==6)
        printtime+=p->P30;
        
        if(p->A10==2 || p->A10==5)
        printtime+=p->P182;
        
        char path[300];
        
        if(p->A10==2)
        sprintf(path,"./REEF3D_SFLOW_6DOF_STL/REEF3D-6DOF-%i-%06i.stl",n6DOF,num);
        
        if(p->A10==5)
        sprintf(path,"./REEF3D_NHFLOW_6DOF_STL/REEF3D-6DOF-%i-%06i.stl",n6DOF,num);
        
        if(p->A10==6)
        sprintf(path,"./REEF3D_CFD_6DOF_STL/REEF3D-6DOF-%i-%06i.stl",n6DOF,num);


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


