/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2026 Tobias Martin

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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include<sys/stat.h>
#include"net_barDyn.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"lagoon_nets.h"
#include<vector>


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
    
	// Print net.vtk
    if(p->mpirank==0 && (((p->count%p->P20==0) && p->P30<0.0)  || (p->simtime>printtime && p->P30>0.0)   || p->count==0))
	{
		printtime += p->P30;

		// the tension at each knot: the mean of its bars'
		std::vector<double> knot_tension(nK, 0.0);
		for (int n = 0; n < nK; ++n)
		{
			double output = 0.0;
			int index = 0;

			for (int i = 0; i < nf; ++i)
			{
				if (Pi[i]==n || Ni[i]==n)
				{
					output += T_(i);
					index++;
				}
			}

			knot_tension[n] = index > 0 ? output/index : 0.0;
		}

		// P 18: the net in the LAGOON store as well (lagoon_nets); P 18 1: instead of the .vtk file
		lagoon_nets::knots k;
		k.xyz.resize(3*nK);
		k.velocity.resize(3*nK);
		k.acceleration.resize(3*nK);
		k.tension.resize(nK);
		for (int n = 0; n < nK; ++n)
		{
			for (int c = 0; c < 3; ++c)
			{
				k.xyz[3*n+c] = float(x_(n,c));
				k.velocity[3*n+c] = float(xdot_(n,c));
				k.acceleration[3*n+c] = float(xdotdot_(n,c));
			}
			k.tension[n] = float(knot_tension[n]);
		}
		for (int i = 0; i < nf; i++)
		{
			k.bars.push_back(Pi[i]);
			k.bars.push_back(Ni[i]);
		}
		for (int i = 0; i < nbK; i++)
		{
			k.bars.push_back(Pb[i]);
			k.bars.push_back(Nb[i]);
		}
		const bool stored = lagoon_nets::output(p, nNet, num, std::move(k));

		if(!(stored && p->P18==1))
		{
		sprintf(name,"./REEF3D_CFD_6DOF_Net/REEF3D-Net-%08i-%06i.vtk",nNet,num);


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
		
		for (int n = 0; n < nK; ++n)
		{
            result<<knot_tension[n]<<endl;
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
}


void net_barDyn::buildNet_cyl(lexer *p)
{
}

void net_barDyn::buildNet_wall(lexer *p)
{     
}
