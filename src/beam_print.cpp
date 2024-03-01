/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2024 Tobias Martin

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
#include"beam.h"
#include"lexer.h"

void beam::print(lexer *p)
{
	int num=0;
	
	if(p->P15==1)
    num = p->printcount_sixdof;

    if(p->P15==2)
    num = p->count;
	
	if(num<0)
	num=0;
    
    // Print beam
	if
	(
		p->mpirank==0 && (((p->count%p->P20==0) && p->P30<0.0)  
		|| (p->simtime>printtime && p->P30>0.0)   
		|| p->count==0)
	)
	{
		printtime+=p->P30;
		
		sprintf(name,"./REEF3D_CFD_Beam/REEF3D-Beam-%i-%06i.vtk",nBeam,num);

		ofstream result;
		result.open(name, ios::binary);
		
		result << "# vtk DataFile Version 2.0" << endl;
		result << "Beam " << nBeam << endl;
		result << "ASCII \nDATASET UNSTRUCTURED_GRID" << endl;
		result << "POINTS " << Ne + 1 << " float" <<endl;
		
		for(int n=0; n<Ne+1; ++n)
		{
			result<<c(0,n)<<" "<<c(1,n)<<" "<<c(2,n)<<endl;
		}
		
		result << "\nCELLS " << Ne << " " << (Ne)*3 <<endl;	
		
		for(int n=0; n<Ne; ++n)
		{
			result<<"2 "<< n << " " << n+1 << endl;
		}
		
		result << "\nCELL_TYPES " << Ne << endl;	
		
		for(int n=0; n<Ne; ++n)
		{
			result<<"3"<<endl;
		}	

		result<<"\nPOINT_DATA " << Ne + 1 <<endl;
		result<<"VECTORS tensLoc double"<<endl;
		for (int n = 0; n < Ne + 1; n++)
		{
			result<<f0(1,n)<<" "<<f0(2,n)<<" "<<f0(3,n)<<endl;
		}
		
        result<<"VECTORS tensGlob double"<<endl;
        Eigen::Vector3d tension;
		for (int n = 0; n < Ne + 1; n++)
		{
            tension = R(q.col(n+1))*f0.col(n).tail(3);
			result<<tension(0)<<" "<<tension(1)<<" "<<tension(2)<<endl;
		}
		
        result<<"VECTORS moments double"<<endl;
		for (int n = 0; n < Ne + 1; n++)
		{
			result<<m0(1,n)<<" "<<m0(2,n)<<" "<<m0(3,n)<<endl;
		}

		result<<"VECTORS velocity double"<<endl;
		for (int n = 0; n < Ne + 1; n++)
		{
			result<<cdot(0,n)<<" "<<cdot(1,n)<<" "<<cdot(2,n)<<endl;
		}
		
        result<<"SCALARS qr float 1 \nLOOKUP_TABLE default"<<endl;
		for (int n = 0; n < Ne + 1; n++)
		{
		    result<<q(0,n)<<endl;
		}
		
		result<<"VECTORS qimg double"<<endl;
		for (int n = 0; n < Ne + 1; n++)
		{
			result<<q(1,n)<<" "<<q(2,n)<<" "<<q(3,n)<<endl;
		}
		
		result<<"VECTORS fext double"<<endl;
		for (int n = 0; n < Ne + 1; n++)
		{
		    result<<Fext(0,n)<<" "<<Fext(1,n)<<" "<<Fext(2,n)<<endl;
		}
		
        result.close();
	}
}

