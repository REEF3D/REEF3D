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
#include"net_barQuasiStatic.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"


void net_barQuasiStatic::print(lexer *p)
{
	int num=0;
	
	if(p->P15==1)
    num = p->printcount_sixdof;

    if(p->P15==2)
    num = p->count;
	
	if(num<0)
	num=0;
	
    // Build current net
    if (p->X320_type[nNet]==1)
    {
        buildNet_bag(p);
    }
    else if (p->X320_type[nNet]==2)
    {
        buildNet_cyl(p);        
    }
    else if (p->X320_type[nNet]==3)
    {
        buildNet_wall(p);        
    }    
    
    // Print maximum tension forces
    if (p->mpirank==0)
    {
        double Tmax = 0.0;
        
		for (int j = 0; j < nf; j++)	
		{
			for (int i = 0; i < niK; i++)	
			{	
				if (fabs(A(i,j)) > Tmax) Tmax = fabs(A(i,j)); 
			}
		}
        
		eTout<<p->simtime<<"\t"<<Tmax<<"\t"<<Fx<<"\t"<<Fy<<"\t"<<Fz<<endl;
    }    
    
	if
	(
		p->mpirank==0 && (((p->count%p->P20==0) && p->P30<0.0)  
		|| (p->simtime>printtime && p->P30>0.0)   
		|| p->count==0)
	)
	{
		printtime+=p->P30;
		
		sprintf(name,"./REEF3D_CFD_6DOF_Net/REEF3D-Net-%08i-%06i.vtk",nNet,num);


		// Tension forces
		double *T;
		p->Darray(T, nf);
	
		for (int j = 0; j < nf; j++)	
		{
			T[j] = 0.0;
			
			for (int i = 0; i < niK; i++)	
			{	
				T[j] = MAX(fabs(A(i,j)),T[j]); 
			}
		}

		ofstream result;
		result.open(name);
		
		result << "# vtk DataFile Version 2.0" << endl;
		result << "Net " << nNet << endl;
		result << "ASCII \nDATASET UNSTRUCTURED_GRID" << endl;
		result << "POINTS " << nK << " float" <<endl;
		
		for(int n=0; n<nK; ++n)
		{
			result<<K_[n][0]<<" "<<K_[n][1]<<" "<<K_[n][2]<<endl;
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
		for (int n = 0; n < nK; ++n)
		{
			output = 0.0;
			int index = 0;
			for (int i = 0; i < nf; ++i)
			{
				if (Pi[i]==n || Ni[i]==n)
				{
					output += T[i]/2.0;
					index++;
				}
			}
			
			result<<output<<endl;
		}
		
		p->del_Darray(T, nf);

		result.close();
	}
}

void net_barQuasiStatic::buildNet_bag(lexer *p)
{
	int *fillK;
	p->Iarray(fillK, nK);
    
    
	for (int i = 0; i < 2*nd+2*nl; i++)
	{
		K_[Pb[i]][0] = K[Pb[i]][0];
		K_[Pb[i]][1] = K[Pb[i]][1];
		K_[Pb[i]][2] = K[Pb[i]][2];			

		K_[Nb[i]][0] = K[Nb[i]][0];
		K_[Nb[i]][1] = K[Nb[i]][1];
		K_[Nb[i]][2] = K[Nb[i]][2];	

		fillK[Pb[i]] = 1;
		fillK[Nb[i]] = 1;
	} 
  
 	K_[nd+1][0] = K[0][0] + fi(0,0)*al;
	K_[nd+1][1] = K[0][1] + fi(0,1)*al;
	K_[nd+1][2] = K[0][2] + fi(0,2)*al;
  
  
	int curRow = 0;
	int leftK = nd + 1;

	for (int i = 0; i < 50; i++)
	{
		if (curRow%2==0)   // even index
		{
			for (int mI = 0; mI < nd - 1; mI++)
			{
				for (int fI = 0; fI < nf; fI++)
				{
					if (Pi[fI]==leftK && Ni[fI]==leftK + nd + 1 && fillK[leftK + nd + 1]==0)
					{
						K_[leftK + nd + 1][0] = K_[leftK][0] + fi(fI,0)*al;	
						K_[leftK + nd + 1][1] = K_[leftK][1] + fi(fI,1)*al;
						K_[leftK + nd + 1][2] = K_[leftK][2] + fi(fI,2)*al;
	
						fillK[leftK + nd + 1] = 1;
					}
					else if (Ni[fI]==leftK && Pi[fI]==leftK + nd + 1 && fillK[leftK + nd + 1]==0)
					{
						K_[leftK + nd + 1][0] = K_[leftK][0] - fi(fI,0)*al;	
						K_[leftK + nd + 1][1] = K_[leftK][1] - fi(fI,1)*al;
						K_[leftK + nd + 1][2] = K_[leftK][2] - fi(fI,2)*al;	

						fillK[leftK + nd + 1] = 1;		
					}
				}

				for (int fI = 0; fI < nf; fI++)
				{
					if (Pi[fI]==leftK + nd + 1 && Ni[fI]==leftK + 1 && fillK[leftK + 1]==0)
					{
						K_[leftK + 1][0] = K_[leftK + nd + 1][0] + fi(fI,0)*al;	
						K_[leftK + 1][1] = K_[leftK + nd + 1][1] + fi(fI,1)*al;
						K_[leftK + 1][2] = K_[leftK + nd + 1][2] + fi(fI,2)*al;

						fillK[leftK + 1] = 1;				
					}
					else if (Ni[fI]==leftK + nd + 1 && Pi[fI]==leftK + 1 && fillK[leftK + 1]==0)
					{
						K_[leftK + 1][0] = K_[leftK + nd + 1][0] - fi(fI,0)*al;	
						K_[leftK + 1][1] = K_[leftK + nd + 1][1] - fi(fI,1)*al;
						K_[leftK + 1][2] = K_[leftK + nd + 1][2] - fi(fI,2)*al;

						fillK[leftK + 1] = 1;				
					}
				}
				leftK++;
			}
			curRow++;
			leftK++;
		}
		else		// odd index
		{
			for (int mI = 0; mI < nd; mI++)
			{
				for (int fI = 0; fI < nf; fI++)
				{
					if (Pi[fI]==leftK && Ni[fI]==leftK + nd + 1 && fillK[leftK + nd + 1]==0)
					{
						K_[leftK + nd + 1][0] = K_[leftK][0] + fi(fI,0)*al;	
						K_[leftK + nd + 1][1] = K_[leftK][1] + fi(fI,1)*al;
						K_[leftK + nd + 1][2] = K_[leftK][2] + fi(fI,2)*al;
	
						fillK[leftK + nd + 1] = 1;
					}
					else if (Ni[fI]==leftK && Pi[fI]==leftK + nd + 1 && fillK[leftK + nd + 1]==0)
					{
						K_[leftK + nd + 1][0] = K_[leftK][0] - fi(fI,0)*al;	
						K_[leftK + nd + 1][1] = K_[leftK][1] - fi(fI,1)*al;
						K_[leftK + nd + 1][2] = K_[leftK][2] - fi(fI,2)*al;	

						fillK[leftK + nd + 1] = 1;		
					}
				}

				for (int fI = 0; fI < nf; fI++)
				{
					if (Pi[fI]==leftK + nd + 1 && Ni[fI]==leftK + 1 && fillK[leftK + 1]==0)
					{
						K_[leftK + 1][0] = K_[leftK + nd + 1][0] + fi(fI,0)*al;	
						K_[leftK + 1][1] = K_[leftK + nd + 1][1] + fi(fI,1)*al;
						K_[leftK + 1][2] = K_[leftK + nd + 1][2] + fi(fI,2)*al;

						fillK[leftK + 1] = 1;				
					}
					else if (Ni[fI]==leftK + nd + 1 && Pi[fI]==leftK + 1 && fillK[leftK + 1]==0)
					{
						K_[leftK + 1][0] = K_[leftK + nd + 1][0] - fi(fI,0)*al;	
						K_[leftK + 1][1] = K_[leftK + nd + 1][1] - fi(fI,1)*al;
						K_[leftK + 1][2] = K_[leftK + nd + 1][2] - fi(fI,2)*al;

						fillK[leftK + 1] = 1;				
					}
				}
				leftK++;
			}
			curRow++;
			leftK++;			
		}
	}
	
	p->del_Iarray(fillK, nK);
}

void net_barQuasiStatic::buildNet_cyl(lexer *p)
{
    // Boundary knots do not move
    for (int i = 0; i < nd; i++)
    {
        K_[i][0] = K[i][0];
        K_[i][1] = K[i][1];
        K_[i][2] = K[i][2];
    }
    
    // Find position of inner knots by going along vertical lines and knowing that owner is above neighbour
    int index = 0;
        
    for (int i = 0; i < nd*nl; i++)
    {
        K_[Ni[i]][0] = K_[Pi[i]][0] + fi(index,0)*l[index];
        K_[Ni[i]][1] = K_[Pi[i]][1] + fi(index,1)*l[index];
        K_[Ni[i]][2] = K_[Pi[i]][2] + fi(index,2)*l[index];
        
        index++;
    }    
}

void net_barQuasiStatic::buildNet_wall(lexer *p)
{
    // Boundary knots do not move
    for (int i = 0; i < nd + 1; i++)
    {
        K_[i][0] = K[i][0];
        K_[i][1] = K[i][1];
        K_[i][2] = K[i][2];
    }
    
    // Find position of inner knots by going along vertical lines and knowing that owner is above neighbour
    int index = 0;
        
    for (int i = 0; i < (nd + 1)*nl; i++)
    {
        K_[Ni[i]][0] = K_[Pi[i]][0] + fi(index,0)*l[index];
        K_[Ni[i]][1] = K_[Pi[i]][1] + fi(index,1)*l[index];
        K_[Ni[i]][2] = K_[Pi[i]][2] + fi(index,2)*l[index];

        index++;
    }  
}
