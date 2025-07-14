/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2025 Tobias Martin

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
Author: Tobias Martin, Hans Bihs
--------------------------------------------------------------------*/

#include"net_barQuasiStatic.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void net_barQuasiStatic::update_velocity_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{

    //- Get velocities at knots
    updateField_nhflow(p, d, pgc, 0);
    updateField_nhflow(p, d, pgc, 1);	
    updateField_nhflow(p, d, pgc, 2);
    
    //- Get density at knots
    updateField_nhflow(p, d, pgc, 3);  
}

void net_barQuasiStatic::updateField_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, int cmp)
{
	int *recField, *count;
	
	p->Iarray(count,p->mpi_size);
	p->Iarray(recField, nK);

	// Get velocities on own processor
	for (int i = 0; i < nK; i++)
	{	
		if 
		(
			K_[i][0] >= xstart[p->mpirank] && K_[i][0] < xend[p->mpirank] &&
			K_[i][1] >= ystart[p->mpirank] && K_[i][1] < yend[p->mpirank] &&
			K_[i][2] >= zstart[p->mpirank] && K_[i][2] < zend[p->mpirank]
		)
		{
			if (cmp==0)
			{
				coupledField[i][cmp] = p->ccipol4V(d->U,d->WL,d->bed,K_[i][0],K_[i][1],K_[i][2]);
			}
			else if (cmp==1)
			{
				coupledField[i][cmp] = p->ccipol4V(d->V,d->WL,d->bed,K_[i][0],K_[i][1],K_[i][2]);
			}
			else if (cmp==2)
			{
				coupledField[i][cmp] =  p->ccipol4V(d->W,d->WL,d->bed,K_[i][0],K_[i][1],K_[i][2]);
			}
			else if (cmp==3)
			{
				coupledField[i][cmp] = p->W1;
			}
            
			recField[i] = -1;
			count[p->mpirank]++;
		}
		else
		{
			for (int j = 0; j < p->mpi_size; j++)
			{	
				if 
				(
					K_[i][0] >= xstart[j] && K_[i][0] < xend[j] &&
					K_[i][1] >= ystart[j] && K_[i][1] < yend[j] &&
					K_[i][2] >= zstart[j] && K_[i][2] < zend[j]
				)
				{
					recField[i] = j;
					count[j]++;
					break;
				}
				else
				{
					recField[i] = -2;
				}
			}			
		}
	}

	
	// Fill array for sending
	double *sendField;
	p->Darray(sendField, count[p->mpirank]);
	
	int counts = 0;
	for (int i = 0; i < nK; i++)
	{
		if (recField[i]==-1)
		{
			sendField[counts] = coupledField[i][cmp];
			counts++;
		}
	}


	// Prepare arrays for receiving
	double **recvField;

	recvField = new double*[p->mpi_size];

	for (int n = 0; n < p->mpi_size; ++n)
	{
		recvField[n] = new double[count[n]];
		
		for (int m = 0; m < count[n]; ++m)
		{
			recvField[n][m] = 0.0;
		}
	}

	
	// Send and receive
	vector<MPI_Request> sreq(p->mpi_size, MPI_REQUEST_NULL);
	vector<MPI_Request> rreq(p->mpi_size, MPI_REQUEST_NULL);
	MPI_Status status;
	
	for (int j = 0; j < p->mpi_size; j++)
	{
		if (j!=p->mpirank)
		{
			if (count[p->mpirank] > 0)
			{
			//	cout<<"Processor "<<p->mpirank<<" sends "<<count[p->mpirank]<<" elements to processor "<<j<<endl;
				
				MPI_Isend(sendField,count[p->mpirank],MPI_DOUBLE,j,1,pgc->mpi_comm,&sreq[j]);
			}
			
			if (count[j] > 0)
			{
			//	cout<<"Processor "<<p->mpirank<<" receives "<<count[j]<<" elements from processor "<<j<<endl;					
		
				MPI_Irecv(recvField[j],count[j],MPI_DOUBLE,j,1,pgc->mpi_comm,&rreq[j]);
			}
		}
	}

	// Wait until transmitted
	for (int j = 0; j < p->mpi_size; j++)
	{
		MPI_Wait(&sreq[j],&status);
		MPI_Wait(&rreq[j],&status);
	}
	
	
	// Fill velocity vector
	for (int j = 0; j < p->mpi_size; j++)
	{
		if (j!=p->mpirank)
		{
			count[j] = 0;
		}
	}
		
	for (int i = 0; i < nK; i++)
	{
		for (int j = 0; j < p->mpi_size; j++)
		{			
			if (recField[i]==j)
			{		
				coupledField[i][cmp] = recvField[j][count[j]];
				count[j]++;
			}
		}
	}
	
	for (int i = 0; i < nK; i++)
	{	 
		coupledField[i][cmp] += 1e-10;
	}	


	// Delete arrays
	if (count[p->mpirank] > 0)
	{
		p->del_Darray(sendField, count[p->mpirank]);
	}

    for(int i = 0; i < p->mpi_size; ++i)
	{
		if (count[i] > 0)
		{
			delete [ ] recvField[i];
		}
	}
	delete [ ] recvField;
	
	p->del_Iarray(count,p->mpi_size);
	p->del_Iarray(recField, nK);
}
