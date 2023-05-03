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

#include"net_sheet.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"reinidisc_fsf.h"	

net_sheet::net_sheet(int number, lexer *p):nNet(number){}

net_sheet::~net_sheet(){}


void net_sheet::initialize(lexer *p, fdm *a, ghostcell *pgc)
{    
    //- Initialise net model
    ini(p,a,pgc);
    
    //- Initialise printing
    printtime = 0.0;
    print(p);
}


void net_sheet::start
(
	lexer *p, 
	fdm *a, 
	ghostcell *pgc,
    double alpha,
    Eigen::Matrix3d quatRotMat
)
{
    double starttime1 = pgc->timer();    
    dt_ = alpha*p->dt;

    //- Store old velocities
    for (int i = 0; i < nK; i++)
    {
        coupledFieldn[i][0] = coupledField[i][0];
        coupledFieldn[i][1] = coupledField[i][1];
        coupledFieldn[i][2] = coupledField[i][2];
    }

    //- Get velocities at knots
    updateField(p, a, pgc, 0);
    updateField(p, a, pgc, 1);	
    updateField(p, a, pgc, 2);
    
    //- Get density at knots
    updateField(p, a, pgc, 3);        
    
    //- Calculate velocities from rigid body motion
    for (int knotI = 0; knotI < nK; knotI++)
    {
        xdot_(knotI,0) = p->ufbi + (x_(knotI,2) - p->zg)*p->qfbi - (x_(knotI,1) - p->yg)*p->rfbi;
        xdot_(knotI,1) = p->vfbi + (x_(knotI,0) - p->xg)*p->rfbi - (x_(knotI,2) - p->zg)*p->pfbi;
        xdot_(knotI,2) = p->wfbi + (x_(knotI,1) - p->yg)*p->pfbi - (x_(knotI,0) - p->xg)*p->qfbi;
    }
    
    //- Calculate force vector
    forces_knot *= 0.0; 
    gravityForce(p);
    inertiaForce(p);
    dragForce(p);
    
    Fx = 0.0;
    Fy = 0.0;
    Fz = 0.0;
    for (int knotI = 0; knotI < nK; knotI++)
    {
        Fx += forces_knot(knotI, 0); 
        Fy += forces_knot(knotI, 1); 
        Fz += forces_knot(knotI, 2); 
    }

	// Update position of triangles
    Eigen::Vector3d point;
	for(int n = 0; n < nK; ++n)
	{
        for(int q = 0; q < 3; q++)
        {
            // (tri_x0 is initial vector between tri_x and xg)
            point << tri_x0[n][q], tri_y0[n][q], tri_z0[n][q];
					
            point = quatRotMat*point;
        
            tri_x[n][q] = point(0) + p->xg;
            tri_y[n][q] = point(1) + p->yg;
            tri_z[n][q] = point(2) + p->zg;
        }
        
        // (x0_ is initial vector between x_ and xg)
        point = quatRotMat*(x0_.row(n)).transpose();
    
        x_(n,0) = point(0) + p->xg;
        x_(n,1) = point(1) + p->yg;
        x_(n,2) = point(2) + p->zg;

	}

    //- Coupling to vrans model
    vransCoupling(p,a,pgc);
	
    //- Build and save net
	print(p);	

    //- Print output
    double endtime1 = pgc->timer() - starttime1; 
    if (p->mpirank == 0)
    {
        cout<<"Net time: "<<endtime1<<endl;    
    }
}


void net_sheet::updateField(lexer *p, fdm *a, ghostcell *pgc, int cmp)
{
	int *recField, *count;

	p->Iarray(count,p->mpi_size);
	p->Iarray(recField, nK);
	
	// Get velocities on own processor
	for (int i = 0; i < nK; i++)
	{	
		if 
		(
			x_(i,0) >= xstart[p->mpirank] && x_(i,0) < xend[p->mpirank] &&
			x_(i,1) >= ystart[p->mpirank] && x_(i,1) < yend[p->mpirank] &&
			x_(i,2) >= zstart[p->mpirank] && x_(i,2) < zend[p->mpirank]
		)
		{
			if (cmp == 0)
			{
				coupledField[i][cmp] = p->ccipol1_a(a->u,x_(i,0),x_(i,1),x_(i,2));
			}
			else if (cmp == 1)
			{
				coupledField[i][cmp] = p->ccipol2_a(a->v,x_(i,0),x_(i,1),x_(i,2));
			}
			else if (cmp == 2)
			{
				coupledField[i][cmp] = p->ccipol3_a(a->w,x_(i,0),x_(i,1),x_(i,2));
			}
			else if (cmp == 3)
			{
				coupledField[i][cmp] = p->ccipol4_a(a->phi,x_(i,0),x_(i,1),x_(i,2));
                
                if (coupledField[i][cmp] >= 0.0) // water
                {
                    coupledField[i][cmp] = p->W1;
                }
                else    // air
                {
                    coupledField[i][cmp] = p->W3;
		}
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
					x_(i,0) >= xstart[j] && x_(i,0) < xend[j] &&
					x_(i,1) >= ystart[j] && x_(i,1) < yend[j] &&
					x_(i,2) >= zstart[j] && x_(i,2) < zend[j]
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
		if (recField[i] == -1)
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
		if (j != p->mpirank)
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
		if (j != p->mpirank)
		{
			count[j] = 0;
		}
	}
		
	for (int i = 0; i < nK; i++)
	{
		for (int j = 0; j < p->mpi_size; j++)
		{			
			if (recField[i] == j)
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


