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

#include"mooring_dynamic.h"
#include"lexer.h"
#include"ghostcell.h"

mooring_dynamic::mooring_dynamic(int number):line(number),beam(number)
{}

mooring_dynamic::~mooring_dynamic(){}

void mooring_dynamic::start(lexer *p, ghostcell *pgc)
{
	// Set mooring time step
	phi_mooring = 0.0;
	t_mooring_n = t_mooring;
	t_mooring = phi_mooring*p->simtime + (1.0 - phi_mooring)*(p->simtime + p->dt);

	// Update fields
    updateFields(p, pgc);
	
	// Update boundary conditions
    fixPoint << p->X311_xe[line], p->X311_ye[line], p->X311_ze[line]; 

    // Integrate from t_mooring_n to t_mooring
    Integrate(t_mooring_n,t_mooring);

	// Save mooring point
	saveMooringPoint(p);

	// Plot mooring line	
	print(p);
}

void mooring_dynamic::updateFluidVel(lexer *p, ghostcell *pgc, int cmp)
{
    /*
	int *recVel, *count;
	
	p->Iarray(count,p->mpi_size);
	p->Iarray(recVel,Ne + 1);
	
	// Get velocities on own processor
	for (int i = 0; i < Ne + 1; i++)
	{	
        fluid_vel[i][cmp] = 0.0;

		if 
		(
			c_moor(0,i) >= xstart[p->mpirank] && c_moor(0,i) < xend[p->mpirank] &&
			c_moor(1,i) >= ystart[p->mpirank] && c_moor(1,i) < yend[p->mpirank] &&
			c_moor(2,i) >= zstart[p->mpirank] && c_moor(2,i) < zend[p->mpirank]
		)
		{
			if (cmp == 0)
			{
				fluid_vel[i][cmp] = p->ccipol1_a(a->u,c_moor(0,i),c_moor(1,i),c_moor(2,i));
			}
			else if (cmp == 1)
			{
				fluid_vel[i][cmp] = p->ccipol2_a(a->v,c_moor(0,i),c_moor(1,i),c_moor(2,i));
			}
			else
			{
				fluid_vel[i][cmp] = p->ccipol3_a(a->w,c_moor(0,i),c_moor(1,i),c_moor(2,i));
			}

			recVel[i] = -1;
			count[p->mpirank]++;
		}
		else
		{
			for (int j = 0; j < p->mpi_size; j++)
			{	
				if 
				(
					c_moor(0,i) >= xstart[j] && c_moor(0,i) < xend[j] &&
					c_moor(1,i) >= ystart[j] && c_moor(1,i) < yend[j] &&
					c_moor(2,i) >= zstart[j] && c_moor(2,i) < zend[j]
				)
				{
					recVel[i] = j;
					count[j]++;
					break;
				}
				else
				{
					recVel[i] = -2;
				}
			}			
		}
	}

	
	// Fill array for sending
	double *sendVel;
	p->Darray(sendVel, count[p->mpirank]);
	
	int counts = 0;
	for (int i = 0; i < Ne + 1; i++)
	{
		if (recVel[i] == -1)
		{
			sendVel[counts] = fluid_vel[i][cmp];
			counts++;
		}
	}


	// Prepare arrays for receiving
	double **recvVel;

	recvVel = new double*[p->mpi_size];

	for (int n = 0; n < p->mpi_size; ++n)
	{
		recvVel[n] = new double[count[n]];
		
		for (int m = 0; m < count[n]; ++m)
		{
			recvVel[n][m] = 0.0;
		}
	}

	
	// Send and receive
	vector<MPI_Request> sreq(p->mpi_size, MPI_REQUEST_NULL);
	vector<MPI_Request> rreq(p->mpi_size, MPI_REQUEST_NULL);
	MPI_Status status;
	
	for (int j = 0; j < p->mpi_size; j++)
	{
		if (j != p->mpirank)
		{
			if (count[p->mpirank] > 0)
			{
			//	cout<<"Processor "<<p->mpirank<<" sends "<<count[p->mpirank]<<" elements to processor "<<j<<endl;
				
				MPI_Isend(sendVel,count[p->mpirank],MPI_DOUBLE,j,1,pgc->mpi_comm,&sreq[j]);
			}
			
			if (count[j] > 0)
			{
			//	cout<<"Processor "<<p->mpirank<<" receives "<<count[j]<<" elements from processor "<<j<<endl;					
		
				MPI_Irecv(recvVel[j],count[j],MPI_DOUBLE,j,1,pgc->mpi_comm,&rreq[j]);
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
		if (j != p->mpirank)
		{
			count[j] = 0;
		}
	}
		
	for (int i = 0; i < Ne + 1; i++)
	{
		for (int j = 0; j < p->mpi_size; j++)
		{			
			if (recVel[i] == j)
			{		
				fluid_vel[i][cmp] = recvVel[j][count[j]];
				count[j]++;
			}
		}
	}
	
	for (int i = 0; i < Ne + 1; i++)
	{	 
		fluid_vel[i][cmp] += 1e-10;		
	}	


	// Delete arrays
	if (count[p->mpirank] > 0)
	{
		p->del_Darray(sendVel, count[p->mpirank]);
	}

    for(int i = 0; i < p->mpi_size; ++i)
	{
		if (count[i] > 0)
		{
			delete [ ] recvVel[i];
		}
	}
	delete [ ] recvVel;
	
	p->del_Iarray(count,p->mpi_size);
	p->del_Iarray(recVel,Ne + 1);*/
}


void mooring_dynamic::updateFields(lexer *p, ghostcell *pgc)
{
    // Get position of points
    getTransPos(c_moor);
	
    // Fluid velocity
	updateFluidVel(p, pgc, 0);
	updateFluidVel(p, pgc, 1);
	updateFluidVel(p, pgc, 2);	
	
    // Fluid acceleration
	for (int i = 0; i < Ne + 1; i++)
	{
        fluid_acc[i][0] = (fluid_vel[i][0] - fluid_vel_n[i][0])/p->dt;
        fluid_acc[i][1] = (fluid_vel[i][1] - fluid_vel_n[i][1])/p->dt;
        fluid_acc[i][2] = (fluid_vel[i][2] - fluid_vel_n[i][2])/p->dt;
	}
}


void mooring_dynamic::mooringForces
(
	double& Xme, double& Yme, double& Zme
)
{
    // Tension forces if line is not broken
    if (broken == false)
    {
        Xme = Xme_; 
        Yme = Yme_;
        Zme = Zme_;
    }

    // Breakage due to max tension force
    if (breakTension > 0.0 && fabs(getTensLoc(Ne)) >= breakTension)
    {
        Xme = 0.0; 
        Yme = 0.0;
        Zme = 0.0;

        broken = true;
    }

    // Breakage due to time limit
    if (breakTime > 0.0 && t_mooring >= breakTime)
    {
        Xme = 0.0; 
        Yme = 0.0;
        Zme = 0.0;

        broken = true;
    }
}


void mooring_dynamic::saveMooringPoint(lexer *p)
{
    // Get position and velocity of points
    getTransPos(c_moor);
    getTransVel(cdot_moor);

    // Save location of line 
    c_moor_n = c_moor;

	// Save acceleration of line 
    cdotdot_moor = (cdot_moor - cdot_moor_n)/p->dt;
    cdot_moor_n = cdot_moor;

	// Save reaction forces at mooring point
    if (p->mpirank==0)
    {
		eTout<<p->simtime<<" \t "<<getTensLoc(Ne)<<endl;
    }
    
    Eigen::Vector3d tension = getTensGlob(Ne);
    Xme_ = -tension(0);
	Yme_ = -tension(1);
	Zme_ = -tension(2);
}
