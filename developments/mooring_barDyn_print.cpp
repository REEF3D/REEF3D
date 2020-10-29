/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2020 Tobias Martin

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
#include"mooring_barDyn.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"


void mooring_barDyn::print(lexer *p)
{
	int num=0;
	
	if(p->P15==1)
    num = p->printcount_sixdof;

    if(p->P15==2)
    num = p->count;
	
	if(num<0)
	num=0;
	
    
    // Print top tension forces
    if (p->mpirank==0)
    {
		eTout<<p->simtime<<" \t "<<T_(nK)<<endl;
    }
    
    // Print line
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
			sprintf(name,"./REEF3D_6DOF_Mooring/REEF3D-Mooring-%d-00000%d.vtk",line,num);

			if(num<100&&num>9)
			sprintf(name,"./REEF3D_6DOF_Mooring/REEF3D-Mooring-%d-0000%d.vtk",line,num);

			if(num<1000&&num>99)
			sprintf(name,"./REEF3D_6DOF_Mooring/REEF3D-Mooring-%d-000%d.vtk",line,num);

			if(num<10000&&num>999)
			sprintf(name,"./REEF3D_6DOF_Mooring/REEF3D-Mooring-%d-00%d.vtk",line,num);

			if(num<100000&&num>9999)
			sprintf(name,"./REEF3D_6DOF_Mooring/REEF3D-Mooring-%d-0%d.vtk",line,num);

			if(num>99999)
			sprintf(name,"./REEF3D_6DOF_Mooring/REEF3D-Mooring-%d-%d.vtk",line,num);
		}	

		ofstream result;
		result.open(name, ios::binary);
		
		result << "# vtk DataFile Version 2.0" << endl;
		result << "Mooring line " << line << endl;
		result << "ASCII \nDATASET UNSTRUCTURED_GRID" << endl;
		result << "POINTS " << nK << " float" <<endl;
		
		for(int n=0; n<nK; ++n)
		{
			result<<x_(n,0)<<" "<<x_(n,1)<<" "<<x_(n,2)<<endl;
		}     
		
		result << "\nCELLS " << nK-1 << " " << (nK-1)*3 <<endl;	
		
		for(int n=0; n<nK-1; ++n)
		{
			result<<"2 "<< n << " " << n+1 << endl;
		}
		
		result << "\nCELL_TYPES " << nK-1 << endl;	
		
		for(int n=0; n<nK-1; ++n)
		{
			result<<"3"<<endl;
		}	

		result<<"\nPOINT_DATA " << nK <<endl;
		result<<"SCALARS Tension float 1 \nLOOKUP_TABLE default"<<endl;
		
		for(int n=0; n<nK; ++n)
		{
			result<<T_(n)<<endl;
		}
		result.close();
	}
}
