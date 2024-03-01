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
#include"net_barDyn.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"


void net_barDyn::print(lexer *p)
{
	int num=0;
	
	if(p->P15==1)
    num = p->printcount_sixdof;

    if(p->P15==2)
    num = p->count;
	
	if(num<0)
	num=0;
 
    // Print tension forces
    if (p->mpirank==0)
    {
        char str[1000];

        sprintf(str,"./REEF3D_CFD_6DOF/REEF3D_6DOF_Net_Forces_%i.dat",nNet);
        ofstream header_out;
        header_out.open(str, std::ofstream::out | std::ofstream::app);
		header_out<<p->simtime<<" \t "<<Tne<<" "<<Fx<<" "<<Fy<<" "<<Fz<<endl;
        header_out.close();
    }  
  
    // Print probe points
    if (p->mpirank==0 && p->X324 > 0)
    {
        for (int pp = 0; pp < p->X324; pp++)
        {
            char str[1000];
            sprintf(str,"./REEF3D_CFD_6DOF/REEF3D_6DOF_Net_%i_Point_Probe_%i.dat",nNet,pp+1);
            ofstream header_out;
            header_out.open(str, std::ofstream::out | std::ofstream::app);
            header_out<<p->simtime<<" \t "<<x_(probeKnot(pp),0)<<" \t "<<x_(probeKnot(pp),1)<<" \t "<<x_(probeKnot(pp),2)<<endl;
            header_out.close();
        }
    }  
    
	
    if
	(
		p->mpirank==0 && (((p->count%p->P20==0) && p->P30<0.0)  
		|| (p->simtime>printtime && p->P30>0.0)   
		|| p->count==0)
	)
	{
		printtime += p->P30;

		sprintf(name,"./REEF3D_CFD_6DOF_Net/REEF3D-Net-%08i-%06i.vtk",nNet,num);


        /*
        char str[1000];
        sprintf(str,"./REEF3D_CFD_6DOF_Net/REEF3D_6DOF_Net_Tension_%i_%i.dat",nNet,num);
        ofstream header_out;
        header_out.open(str, std::ofstream::out | std::ofstream::app);
		
        for (int j = 0; j < nf; j++)	
        {    
            header_out<<T_(j)<<" "<<(0.5*(x_.row(Pi[j])+x_.row(Ni[j])))<<endl;
        }
        
        header_out.close();
        */


		ofstream result;
		result.open(name, ios::binary);
		
        result << "# vtk DataFile Version 2.0" << endl;
		result << "Net " << nNet << endl;
		result << "ASCII \nDATASET UNSTRUCTURED_GRID" << endl;
		result << "POINTS " << nK << " float" <<endl;
		
		for(int n=0; n<nK; ++n)
		{
			result<<x_(n,0)<<" "<<x_(n,1)<<" "<<x_(n,2)<<endl;
		}     
        
		result << "\nCELLS " << nf+nbK << " " << (nf+nbK)*3 <<endl;
		
		for (int i = 0; i < nf; i++)
		{
			result<<"2 "<< Pi[i] << " " << Ni[i] <<endl;
		}

		for (int i = 0; i < nbK; i++)
		{
			result<<"2 "<< Pb[i] << " " << Nb[i] <<endl;
		}

		result << "\nCELL_TYPES " << nf+nbK << endl;	
		
		for (int i = 0; i < nf+nbK; i++)
		{
			result<<"3"<<endl;
		}		

		result<<"\nPOINT_DATA " << nK <<endl;
		result<<"SCALARS Tension float 1 \nLOOKUP_TABLE default"<<endl;
		
		double output;
        int index;
        
		for (int n = 0; n < nK; ++n)
		{
			output = 0.0;
            index = 0;
            
			for (int i = 0; i < nf; ++i)
			{
				if (Pi[i] == n || Ni[i] == n)
				{
					output += T_(i);
					index++;
				}
			}
            
            output = index > 0 ? output/index : 0.0;
            result<<output<<endl;
		}

		result<<"VECTORS Velocity float"<<endl;
		for (int i = 0; i < nK; i++)
		{
		    result<<xdot_.row(i)<<endl;
		}

		result<<"VECTORS Acceleration float"<<endl;
		for (int i = 0; i < nK; i++)
		{
		    result<<xdotdot_.row(i)<<endl;
		}

		result.close();
	}
}


void net_barDyn::buildNet_cyl(lexer *p)
{
}

void net_barDyn::buildNet_wall(lexer *p)
{     
}
