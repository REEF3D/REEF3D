/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
#include"mooring_DGSEM.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"



void mooring_DGSEM::print
(
	lexer *p
)
{
	int num=0;
	
	if(p->P15==1)
    num = p->printcount_sixdof;

    if(p->P15==2)
    num = p->count;
	
	if(num<0)
	num=0;
	
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
		result << "POINTS " << H*(P+1) << " double" <<endl;
		
		for (int i = 0; i < H; i++)
		{
			for (int j = 0; j < (P+1); j++)
			{
				result<<r_x[i][j]<<" "<<r_y[i][j]<<" "<<r_z[i][j]<<endl;
			}
		}
		
		result << "\nCELLS " << H << " " << H*3 <<endl;	
		
		for(int n = 0; n < H; ++n)
		{
			result<<"2 "<< n*(P+1) << " " << (n+1)*(P+1)-1 << endl;
		}
		
		result << "\nCELL_TYPES " << H << endl;	
		
		for(int n = 0; n < H; ++n)
		{
			result<<"3"<<endl;
		}	

		result<<"\nCELL_DATA " << H <<endl;		
		result<<"SCALARS SmoothInd double 1"<<endl;
		result<<"LOOKUP_TABLE default"<<endl;
		
		result<<MAX(0,fabs(T[0][P]-T[1][0]))*pow((vx[1]-vx[0]), -0.5*(P+1))<<endl;
		for (int i = 1; i < H-1; i++)
		{
			result<<MAX(fabs(T[i][0]-T[i-1][P]),fabs(T[i][P]-T[i+1][0]))*pow((vx[1]-vx[0]), -0.5*(P+1))<<endl;
		}
		result<<MAX(0,fabs(T[H-1][0]-T[H-2][P]))*pow((vx[1]-vx[0]), -0.5*(P+1))<<endl;


		result<<"SCALARS discError double 1"<<endl;
		result<<"LOOKUP_TABLE default"<<endl;
		
		result<<1.0/sqrt(8.0)*sqrt(pow(fabs(T[0][P]-T[1][0]),2))<<endl;
		for (int i = 1; i < H-1; i++)
		{
			result<<1.0/sqrt(8.0)*sqrt(pow(fabs(T[i][0]-T[i-1][P]),2) + pow(fabs(T[i][P]-T[i+1][0]),2))<<endl;
		}
		result<<1.0/sqrt(8.0)*sqrt(pow(fabs(T[H-1][0]-T[H-2][P]),2))<<endl;
		
	
		result<<"\nPOINT_DATA " << H*(P+1) <<endl;
		
		result<<"VECTORS tension double"<<endl;
		for (int i = 0; i < H; i++)
		{
			for (int j = 0; j < (P+1); j++)
			{
				result<<T[i][j]*t_x[i][j]<<" "<<T[i][j]*t_y[i][j]<<" "<<T[i][j]*t_z[i][j]<<endl;
			}
		}

		result<<"VECTORS velocity double"<<endl;
		for (int i = 0; i < H; i++)
		{
			for (int j = 0; j < (P+1); j++)
			{
				result<<v_x[i][j]<<" "<<v_y[i][j]<<" "<<v_z[i][j]<<endl;
			}
		}
		
		result<<"VECTORS curvature double"<<endl;
		for (int i = 0; i < H; i++)
		{
			for (int j = 0; j < (P+1); j++)
			{
				result<<q_x[i][j]<<" "<<q_y[i][j]<<" "<<q_z[i][j]<<endl;
			}
		}
		
		result<<"VECTORS forces double"<<endl;
		for (int i = 0; i < H; i++)
		{
			for (int j = 0; j < (P+1); j++)
			{
				result<<Fx[i][j]<<" "<<Fy[i][j]<<" "<<Fz[i][j]<<endl;
			}
		}

		result.close();

		eTout<<p->simtime<<" \t "<<T[H-1][P]<<endl;
	}
}