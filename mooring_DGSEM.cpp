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

#include"mooring_DGSEM.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

mooring_DGSEM::mooring_DGSEM(int number):line(number){}

mooring_DGSEM::~mooring_DGSEM(){}

void mooring_DGSEM::start(lexer *p, fdm *a, ghostcell *pgc)
{
	// Set mooring time step
	double phi = 0.5;
	t_mooring_n = t_mooring;
	t_mooring = phi*p->simtime + (1.0 - phi)*(p->simtime + p->dt);
	dtm = t_mooring - t_mooring_n;
	
	dt = cfl*relFac;
	dtau = 0.0;
	
	// Start loop
	if (dt >= dtm)  
	{
		dt = dtm;

		startRK3TVD(p,a,pgc);
	}
	else
	{
		int loops = ceil(dtm/dt);
		dt = dtm/loops;

		for (int loop = 0; loop < loops; loop++)
		{
			startRK3TVD(p,a,pgc);
		}
	}
	
	// Save mooring point
	saveMooringPoint(p);

	// Plot mooring line	
	print(p);
	
	
if (p->mpirank == 0)
{	
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < (P+1); j++)
		{
//			cout<<r_x[i][j]<<" "<<r_z[i][j]<<" "<<T[i][j]<<endl;
		}
	}	
}

}


void mooring_DGSEM::startRK3TVD(lexer *p, fdm *a, ghostcell *pgc)
{
	// Update fluid fields
	updateFields(p, a, pgc);
	
	// Update boundary conditions
	updateBC(p);
	
	
	// Step 1
	getL(p,r_x,r_y,r_z,q_x,q_y,q_z,v_x,v_y,v_z);
	
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < (P+1); j++)
		{
			r_x1[i][j] = r_x[i][j] + dt*Lr_x[i][j];                            
			r_y1[i][j] = r_y[i][j] + dt*Lr_y[i][j];
			r_z1[i][j] = r_z[i][j] + dt*Lr_z[i][j];
			q_x1[i][j] = q_x[i][j] + dt*Lq_x[i][j];
			q_y1[i][j] = q_y[i][j] + dt*Lq_y[i][j];
			q_z1[i][j] = q_z[i][j] + dt*Lq_z[i][j];
			v_x1[i][j] = v_x[i][j] + dt*Lv_x[i][j];
			v_y1[i][j] = v_y[i][j] + dt*Lv_y[i][j];
			v_z1[i][j] = v_z[i][j] + dt*Lv_z[i][j];
		}
	}

	// Apply minMod slope limiter
	limitFields(p,r_x1,r_y1,r_z1,q_x1,q_y1,q_z1,v_x1,v_y1,v_z1);
	
	
	// Step 2
	getL(p,r_x1,r_y1,r_z1,q_x1,q_y1,q_z1,v_x1,v_y1,v_z1);
	
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < (P+1); j++)
		{
			r_x2[i][j] = (3.0*r_x[i][j] + r_x1[i][j] + dt*Lr_x[i][j])/4.0;                            
			r_y2[i][j] = (3.0*r_y[i][j] + r_y1[i][j] + dt*Lr_y[i][j])/4.0;
			r_z2[i][j] = (3.0*r_z[i][j] + r_z1[i][j] + dt*Lr_z[i][j])/4.0;
			q_x2[i][j] = (3.0*q_x[i][j] + q_x1[i][j] + dt*Lq_x[i][j])/4.0;
			q_y2[i][j] = (3.0*q_y[i][j] + q_y1[i][j] + dt*Lq_y[i][j])/4.0;
			q_z2[i][j] = (3.0*q_z[i][j] + q_z1[i][j] + dt*Lq_z[i][j])/4.0;
			v_x2[i][j] = (3.0*v_x[i][j] + v_x1[i][j] + dt*Lv_x[i][j])/4.0;
			v_y2[i][j] = (3.0*v_y[i][j] + v_y1[i][j] + dt*Lv_y[i][j])/4.0;
			v_z2[i][j] = (3.0*v_z[i][j] + v_z1[i][j] + dt*Lv_z[i][j])/4.0;
		}
	}

	// Apply minMod slope limiter
	limitFields(p,r_x2,r_y2,r_z2,q_x2,q_y2,q_z2,v_x2,v_y2,v_z2);
	
	
	// Step 3
	getL(p,r_x2,r_y2,r_z2,q_x2,q_y2,q_z2,v_x2,v_y2,v_z2);
	
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < (P+1); j++)
		{
			r_x[i][j] = (r_x[i][j] + 2.0*r_x2[i][j] + 2.0*dt*Lr_x[i][j])/3.0;                            
			r_y[i][j] = (r_y[i][j] + 2.0*r_y2[i][j] + 2.0*dt*Lr_y[i][j])/3.0;
			r_z[i][j] = (r_z[i][j] + 2.0*r_z2[i][j] + 2.0*dt*Lr_z[i][j])/3.0;
			q_x[i][j] = (q_x[i][j] + 2.0*q_x2[i][j] + 2.0*dt*Lq_x[i][j])/3.0;
			q_y[i][j] = (q_y[i][j] + 2.0*q_y2[i][j] + 2.0*dt*Lq_y[i][j])/3.0;
			q_z[i][j] = (q_z[i][j] + 2.0*q_z2[i][j] + 2.0*dt*Lq_z[i][j])/3.0;
			v_x[i][j] = (v_x[i][j] + 2.0*v_x2[i][j] + 2.0*dt*Lv_x[i][j])/3.0;
			v_y[i][j] = (v_y[i][j] + 2.0*v_y2[i][j] + 2.0*dt*Lv_y[i][j])/3.0;
			v_z[i][j] = (v_z[i][j] + 2.0*v_z2[i][j] + 2.0*dt*Lv_z[i][j])/3.0;
		}
	}
	
	// Apply minMod slope limiter
	limitFields(p,r_x,r_y,r_z,q_x,q_y,q_z,v_x,v_y,v_z);
}


void mooring_DGSEM::updateFluidVel(lexer *p, fdm *a, ghostcell *pgc, int cmp)
{
	int **recVel, *count;
	
	p->Iarray(count,p->mpi_size);
	p->Iarray(recVel, H, P + 1);
	
	// Get velocities on own processor
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < (P+1); j++)
		{
			// Save old fluid velocities
			vfn[i][j][cmp] = vf[i][j][cmp];
			
			if 
			(
				r_x[i][j] >= xstart[p->mpirank] && r_x[i][j] < xend[p->mpirank] &&
				r_y[i][j] >= ystart[p->mpirank] && r_y[i][j] < yend[p->mpirank] &&
				r_z[i][j] >= zstart[p->mpirank] && r_z[i][j] < zend[p->mpirank]
			)
			{
				if (cmp == 0)
				{
					vf[i][j][cmp] = p->ccipol1_a(a->u,r_x[i][j],r_y[i][j],r_z[i][j]);
				}
				else if (cmp == 1)
				{
					vf[i][j][cmp] = p->ccipol2_a(a->v,r_x[i][j],r_y[i][j],r_z[i][j]);
				}
				else
				{
					vf[i][j][cmp] = p->ccipol3_a(a->w,r_x[i][j],r_y[i][j],r_z[i][j]);
				}

				recVel[i][j] = -1;
				count[p->mpirank]++;
			}
			else
			{
				for (int k = 0; k < p->mpi_size; k++)
				{	
					if 
					(
						r_x[i][j] >= xstart[k] && r_x[i][j] < xend[k] &&
						r_y[i][j] >= ystart[k] && r_y[i][j] < yend[k] &&
						r_z[i][j] >= zstart[k] && r_z[i][j] < zend[k]
					)
					{
						recVel[i][j] = k;
						count[k]++;
						break;
					}
					else
					{
						recVel[i][j] = -2;
					}
				}			
			}
		}
	}
	
	// Fill array for sending
	double *sendVel;
	p->Darray(sendVel, count[p->mpirank]);
	
	int counts = 0;
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < (P+1); j++)
		{
			if (recVel[i][j] == -1)
			{
				sendVel[counts] = vf[i][j][cmp];
				counts++;
			}
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
		
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < (P+1); j++)
		{
			for (int k = 0; k < p->mpi_size; k++)
			{			
				if (recVel[i][j] == k)
				{		
					vf[i][j][cmp] = recvVel[k][count[k]];
					count[k]++;
				}
			}
		}
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
	p->del_Iarray(recVel,H, P + 1);
}


void mooring_DGSEM::updateFields(lexer *p, fdm *a, ghostcell *pgc)
{
	// Fluid velocity
	updateFluidVel(p, a, pgc, 0);
	updateFluidVel(p, a, pgc, 1);
	updateFluidVel(p, a, pgc, 2);
	
	
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < (P+1); j++)
		{
			// Fluid acceleration
			af[i][j][0] = (vf[i][j][0] - vfn[i][j][0])/dt;
			af[i][j][1] = (vf[i][j][1] - vfn[i][j][1])/dt;
			af[i][j][2] = (vf[i][j][2] - vfn[i][j][2])/dt;
			
			// Save old mooring velocities
			v_xn[i][j] = v_x[i][j];
			v_yn[i][j] = v_y[i][j];
			v_zn[i][j] = v_z[i][j];
		}
	}
}


void mooring_DGSEM::mooringForces
(
	double& Xme, double& Yme, double& Zme
)
{
	Xme = Xme_; 
	Yme = Yme_;
	Zme = Zme_;
}


void mooring_DGSEM::saveMooringPoint(lexer *p)
{
	// Save location of mooring point
	r_xOn = r_xO;
	r_yOn = r_yO;
	r_zOn = r_zO;

	// Save velocities at mooring point
	v_xOn = v_xO;
	v_yOn = v_yO;
	v_zOn = v_zO;

	// Save reaction forces at mooring point
	Xme_ = -T[H-1][P]*t_x[H-1][P];
	Yme_ = -T[H-1][P]*t_y[H-1][P];
	Zme_ = -T[H-1][P]*t_z[H-1][P];
}
