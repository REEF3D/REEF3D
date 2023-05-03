/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2023 Tobias Martin

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
#include"mooring_Spring.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"


void mooring_Spring::print(lexer *p)
{
	int num=0;
	
	if(p->P15==1)
    num = p->printcount_sixdof;

    if(p->P15==2)
    num = p->count;
	
	if(num<0)
	num=0;
	
    if (p->mpirank==0)
    {
		eTout<<p->simtime<<" \t "<<T<<endl;	
    }  
	
    if
	(
		p->mpirank==0 && (((p->count%p->P20==0) && p->P30<0.0)  
		|| (p->simtime>printtime && p->P30>0.0)   
		|| p->count==0)
	)
	{
		printtime+=p->P30;
		
		if(p->P14==1)
		{
			if(num<10)
			sprintf(name,"./REEF3D_CFD_6DOF_Mooring/REEF3D-Mooring-%i-00000%i.vtk",line,num);

			if(num<100&&num>9)
			sprintf(name,"./REEF3D_CFD_6DOF_Mooring/REEF3D-Mooring-%i-0000%i.vtk",line,num);

			if(num<1000&&num>99)
			sprintf(name,"./REEF3D_CFD_6DOF_Mooring/REEF3D-Mooring-%i-000%i.vtk",line,num);

			if(num<10000&&num>999)
			sprintf(name,"./REEF3D_CFD_6DOF_Mooring/REEF3D-Mooring-%i-00%i.vtk",line,num);

			if(num<100000&&num>9999)
			sprintf(name,"./REEF3D_CFD_6DOF_Mooring/REEF3D-Mooring-%i-0%i.vtk",line,num);

			if(num>99999)
			sprintf(name,"./REEF3D_CFD_6DOF_Mooring/REEF3D-Mooring-%i-%i.vtk",line,num);
		}

		// Print results
		ofstream result;
		result.open(name, ios::binary);
		
		result << "# vtk DataFile Version 2.0" << endl;
		result << "Mooring line " << line << endl;
		result << "ASCII \nDATASET UNSTRUCTURED_GRID" << endl;
		result << "POINTS " << 2 << " float" <<endl;

        result<<p->X311_xs[line]<<" "<<p->X311_ys[line]<<" "<<p->X311_zs[line]<<endl;
        result<<p->X311_xe[line]<<" "<<p->X311_ye[line]<<" "<<p->X311_ze[line]<<endl;
        
		
		result << "\nCELLS " << 1 << " " << 3 <<endl;	
        
        result<<"2 "<< 0 << " " << 1 << endl;

		
		result << "\nCELL_TYPES " << 1 << endl;	

        result<<"3"<<endl;


		result<<"\nPOINT_DATA " << 2 <<endl;
		result<<"SCALARS Tension float 1 \nLOOKUP_TABLE default"<<endl;
		
        result<<T<<endl;
        result<<T<<endl;
		
        
		result.close();


	}
}
