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

#include"mooring_barQuasiStatic.h"
#include"lexer.h"
#include"ghostcell.h"

mooring_barQuasiStatic::mooring_barQuasiStatic(int number):line(number){}

mooring_barQuasiStatic::~mooring_barQuasiStatic(){}

void mooring_barQuasiStatic::start(lexer *p, ghostcell *pgc)
{
    // Current time
    curr_time = p->simtime;

	// Correct geometrical constraint
	dx = p->X311_xe[line] - p->X311_xs[line];
	dy = p->X311_ye[line] - p->X311_ys[line];
	dz = p->X311_ze[line] - p->X311_zs[line]; 
	
	B[sigma][0] = dx;
	B[sigma][1] = dy;
	B[sigma][2] = dz;		

	double length = 0.0;

	// Solving the system of equations
    if (broken==false)
    {
        for (int it = 0; it < 1000; it++)
        {	
            // Reconstruct current line
            buildLine(p,pgc);
            
            // Calculating velocities at knots
            // updateVel(p, a, pgc, 0);
            // updateVel(p, a, pgc, 1);
            // updateVel(p, a, pgc, 2);

            // Calculating force directions
            for (int j = 0; j < sigma + 1; j++)
            {
                v[j][0] = 1e-10;
                v[j][1] = 1e-10;
                v[j][2] = 1e-10;
                
                
                double vt = v[j][0]*f[j][0] + v[j][1]*f[j][1] + v[j][2]*f[j][2];
                double vt_mag = fabs(vt);

                e_d[j][0] = vt_mag*vt*f[j][0];
                e_d[j][1] = vt_mag*vt*f[j][1];
                e_d[j][2] = vt_mag*vt*f[j][2];

                double vvtt_x = v[j][0] - vt*f[j][0];
                double vvtt_y = v[j][1] - vt*f[j][1];
                double vvtt_z = v[j][2] - vt*f[j][2];
                
                vt_mag = sqrt(vvtt_x*vvtt_x + vvtt_y*vvtt_y + vvtt_z*vvtt_z);
                
                e_l[j][0] = vt_mag*vvtt_x;
                e_l[j][1] = vt_mag*vvtt_y;
                e_l[j][2] = vt_mag*vvtt_z;

            /*
                double magVal = sqrt(v[j][0]*v[j][0]+v[j][1]*v[j][1]+v[j][2]*v[j][2]);

                e[j][0] = v[j][0]/magVal;	
                e[j][1] = v[j][1]/magVal;	
                e[j][2] = v[j][2]/magVal;

                // Drag force
                e_d[j][0] = -e[j][0];
                e_d[j][1] = -e[j][1];
                e_d[j][2] = -e[j][2];			
                
                // Shear force
                double val_x = f[j][2]*e[j][1] - f[j][1]*e[j][2];
                double val_y = f[j][0]*e[j][2] - f[j][2]*e[j][0];
                double val_z = f[j][1]*e[j][0] - f[j][0]*e[j][1];
                    
                magVal = sqrt(val_x*val_x+val_y*val_y+val_z*val_z);
                    
                e_q[j][0] = val_x / magVal;
                e_q[j][1] = val_y / magVal;
                e_q[j][2] = val_z / magVal;
                
                //	Lift force	
                val_x = e_q[j][2]*e[j][1] - e_q[j][1]*e[j][2];
                val_y = e_q[j][0]*e[j][2] - e_q[j][2]*e[j][0];
                val_z = e_q[j][1]*e[j][0] - e_q[j][0]*e[j][1];			

                magVal = sqrt(val_x*val_x+val_y*val_y+val_z*val_z);
                    
                e_l[j][0] = val_x / magVal;
                e_l[j][1] = val_y / magVal;
                e_l[j][2] = val_z / magVal;

                // Calculating force coefficients
                if (sqrt(f[j][0]*f[j][0] + f[j][1]*f[j][1] + f[j][2]*f[j][2]) <= 1.0)
                {
                    double theta = acos(e[j][0]*f[j][0] + e[j][1]*f[j][1] + e[j][2]*f[j][2]);
                        
                    c_coeff[j] = getC(theta);
                }	*/		
            }
                
            // Calculating hydrodynamic forces
            for (int j=1; j<sigma+1; j++)
            {
                double cdt = 0.5;
                double cdn = 2.5;

                R[j][0] = 
                    p->W1/2.0*p->X311_d[line]*l[j-1]/2.0*(cdt*e_d[j-1][0] + cdn*e_l[j-1][0])
                    + p->W1/2.0*p->X311_d[line]*l[j]/2.0*(cdt*e_d[j][0] + cdn*e_l[j][0]);
                
                R[j][1] = 
                    p->W1/2.0*p->X311_d[line]*l[j-1]/2.0*(cdt*e_d[j-1][1] + cdn*e_l[j-1][1])
                    + p->W1/2.0*p->X311_d[line]*l[j]/2.0*(cdt*e_d[j][1] + cdn*e_l[j][1]);
                        
                R[j][2] = 
                    p->W1/2.0*p->X311_d[line]*l[j-1]/2.0*(cdt*e_d[j-1][2] + cdn*e_l[j-1][2])
                    + p->W1/2.0*p->X311_d[line]*l[j]/2.0*(cdt*e_d[j][2] + cdn*e_l[j][2]);
                    
                
            /*	R[j][0] = 
                    0.5*p->W1*v[j][0]*v[j][0]*0.5*l[j]*p->X311_d[line]*
                    (
                          c_coeff[j-1][0]*e_d[j-1][0] + c_coeff[j][0]*e_d[j][0] 
                        + c_coeff[j-1][1]*e_q[j-1][0] + c_coeff[j][1]*e_q[j][0] 
                        + c_coeff[j-1][2]*e_l[j-1][0] + c_coeff[j][2]*e_l[j][0]
                    );
                R[j][1] = 
                    0.5*p->W1*v[j][1]*v[j][1]*0.5*l[j]*p->X311_d[line]*
                    (
                          c_coeff[j-1][0]*e_d[j-1][1] + c_coeff[j][0]*e_d[j][1] 
                        + c_coeff[j-1][1]*e_q[j-1][1] + c_coeff[j][1]*e_q[j][1] 
                        + c_coeff[j-1][2]*e_l[j-1][1] + c_coeff[j][2]*e_l[j][1]
                    );	
                R[j][2] = 
                    0.5*p->W1*v[j][2]*v[j][2]*0.5*l[j]*p->X311_d[line]*
                    (
                          c_coeff[j-1][0]*e_d[j-1][2] + c_coeff[j][0]*e_d[j][2] 
                        + c_coeff[j-1][1]*e_q[j-1][2] + c_coeff[j][1]*e_q[j][2] 
                        + c_coeff[j-1][2]*e_l[j-1][2] + c_coeff[j][2]*e_l[j][2]
                    );	*/		
            } 
                
            // Filling right hand side
            for (int j=0; j<sigma; j++)
            {
                B[j][0] = -R[j+1][0];
                B[j][1] = -R[j+1][1];
                B[j][2] = -R[j+1][2] + l[j+1] * w - Fb[j+1];
            }
            
            // Solving the system
            f = solveGauss(A, B);

            // Check error norm
            double norm, error;
                
            error = 0.0;
            for (int j=0; j<sigma+1; j++)
            {
                norm = sqrt(f[j][0]*f[j][0] + f[j][1]*f[j][1] + f[j][2]*f[j][2]);
                
                if (norm > error) error = norm;
            }

            if (fabs(error-1.0) < 1e-4)
            {
                if (p->mpirank==0)
                {
                    cout<<"Number of iterations = "<<it<<setprecision(6)<<" with error = "<<error-1.0<<endl;
                }
                break;
            }
                
            // Correct system
            for (int j = 0; j < sigma + 1; j++)	
            {
                norm = sqrt(f[j][0]*f[j][0] + f[j][1]*f[j][1] + f[j][2]*f[j][2]);
                    
                f[j][0] /= norm;
                f[j][1] /= norm;
                f[j][2] /= norm;
            
                for (int k=0; k<sigma; k++) 
                {
                    A[k][j] *= norm;
                }
            }

            length = 0.0;
            for (int j = 1; j < sigma+1; j++)	
            {			
                l[j] = l0[j]*(1.0 + 0.5*(fabs(A[j][j]) + fabs(A[j-1][j-1]))/(p->X311_EA[line]));
                length += l[j];
            }
                
            for (int j = 0; j < sigma; j++)
            {			
                A[sigma][j] = 0.5*(l[j] + l[j+1]);  
            }
            A[sigma][sigma] = 0.5*l[sigma];
        }	
    }
    else
    {
        if (p->mpirank==0) cout<<"Line "<<line<<" broken"<<endl;
    }

    //if (p->mpirank==0) cout<<"Current length = "<<length<<endl;
	
	// Calculate horizontal and vertical reaction forces at mooring point
	Xme_ = -fabs(A[sigma-1][sigma])*f[sigma][0];
	Yme_ = -fabs(A[sigma-1][sigma])*f[sigma][1];
	Zme_ = -fabs(A[sigma-1][sigma])*f[sigma][2];
		
	// Plotting mooring line	
	print(p,pgc);	
}


void mooring_barQuasiStatic::updateVel(lexer *p, ghostcell *pgc, int cmp)
{
	/*int *recVel, *count;
	
	p->Iarray(count,p->mpi_size);
	p->Iarray(recVel,sigma + 2);
	
	// Get velocities on own processor
	for (int i = 0; i < sigma + 2; i++)
	{	
		if 
		(
			x[i] >= xstart[p->mpirank] && x[i] < xend[p->mpirank] &&
			y[i] >= ystart[p->mpirank] && y[i] < yend[p->mpirank] &&
			z[i] >= zstart[p->mpirank] && z[i] < zend[p->mpirank]
		)
		{
			if (cmp==0)
			{
				v[i][cmp] = p->ccipol1_a(a->u,x[i],y[i],z[i]);
			}
			else if (cmp==1)
			{
				v[i][cmp] = p->ccipol2_a(a->v,x[i],y[i],z[i]);
			}
			else
			{
				v[i][cmp] = p->ccipol3_a(a->w,x[i],y[i],z[i]);
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
					x[i] >= xstart[j] && x[i] < xend[j] &&
					y[i] >= ystart[j] && y[i] < yend[j] &&
					z[i] >= zstart[j] && z[i] < zend[j]
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
	for (int i = 0; i < sigma + 2; i++)
	{
		if (recVel[i]==-1)
		{
			sendVel[counts] = v[i][cmp];
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
		if (j!=p->mpirank)
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
		if (j!=p->mpirank)
		{
			count[j] = 0;
		}
	}
		
	for (int i = 0; i < sigma + 2; i++)
	{
		for (int j = 0; j < p->mpi_size; j++)
		{			
			if (recVel[i]==j)
			{		
				v[i][cmp] = recvVel[j][count[j]];
				count[j]++;
			}
		}
	}
	
	for (int i = 0; i < sigma + 2; i++)
	{	 
		v[i][cmp] += 1e-10;		
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
	p->del_Iarray(recVel,sigma + 2);*/
}


vector<double> mooring_barQuasiStatic::getC(double theta)
{
    vector<double> c_(3,0);
    
	theta *= 180/PI;
	
	c_[0] = max(0.0, -0.000000016676062*pow(theta,4.0) + 0.000001614232301*pow(theta,3.0) + 0.000132218750548*pow(theta,2.0) - 0.001335313221272*theta + 0.023082880162714);
	c_[1] = max(0.0, 0.000000019889794877*pow(theta,4.0) - 0.000005200545357615*pow(theta,3.0) + 0.000309491763289986*pow(theta,2.0) + 0.000866335727454720*theta - 0.000960622162325914);
	c_[2] = max(0.0, -0.000000028334918*pow(theta,4.0) + 0.000000348243732*pow(theta,3.0) + 0.000255407378882*pow(theta,2.0) - 0.006366497021954*theta + 0.049082170481768);

    return c_;
}



vector< vector<double> > mooring_barQuasiStatic::solveGauss
(
	std::vector< std::vector<double> > A, 
	std::vector< std::vector<double> > B
) 
{
    int n  = A.size();
	int nB = B[0].size();

    for (int i=0; i<n; i++) 
	{
        // Search for maximum in this column
        double maxEl = fabs(A[i][i]);
        int maxRow = i;
		
        for (int k=i+1; k<n; k++) 
		{
            if (fabs(A[k][i]) > maxEl) 
			{
                maxEl = fabs(A[k][i]);
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        for (int k=i; k<n; k++)
		{
            double tmp = A[maxRow][k];
			
            A[maxRow][k] = A[i][k];
            A[i][k] = tmp;
        }
        for (int k=0; k<nB; k++)
		{
            double tmp = B[maxRow][k];
			
            B[maxRow][k] = B[i][k];
            B[i][k] = tmp;
        }
		
		
        // Make all rows below this one 0 in current column
        for (int k=i+1; k<n; k++) 
		{
            double c = -A[k][i]/A[i][i];
			
            for (int j=i; j<n; j++)
			{
                if (i==j) 
				{
                    A[k][j] = 0;
				} 
				else 
				{
                    A[k][j] += c * A[i][j];
                }
            }
            for (int j=0; j<nB; j++)
			{
				B[k][j] += c * B[i][j];
            }			
        }
    }

    // Solve equation Ax=b for an upper triangular matrix A
	vector<double> xVec(nB);
    vector< vector<double> > x(n,xVec);

    for (int i=n-1; i>=0; i--) 
	{
		for (int j=0; j<nB; j++)
		{
			x[i][j] = B[i][j]/A[i][i];
		}			
		
        for (int k=i-1; k>=0; k--) 
		{
			for (int j=0; j<nB; j++)
			{
				B[k][j] -= A[k][i] * x[i][j];
			}			
        }
    }

    return x;
}

void mooring_barQuasiStatic::mooringForces
(
	double& Xme, double& Yme, double& Zme
)
{
    // Tension forces if line is not broken
    if (broken==false)
    {
        Xme = Xme_; 
        Yme = Yme_;
        Zme = Zme_;
    }

    // Breakage due to max tension force
    if (breakTension > 0.0 && fabs(A[sigma-1][sigma]) >= breakTension)
    {
        Xme = 0.0; 
        Yme = 0.0;
        Zme = 0.0;

        broken = true;
    }

    // Breakage due to time limit
    if (breakTime > 0.0 && curr_time >= breakTime)
    {
        Xme = 0.0; 
        Yme = 0.0;
        Zme = 0.0;

        broken = true;
    }
}
