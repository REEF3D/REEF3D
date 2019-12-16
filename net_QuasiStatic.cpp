/*--------------------------------------------------------------------
REEF3D
Copyright 2019 Tobias Martin

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

#include"net_QuasiStatic.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
//#include <Eigen/Dense>

//typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;	

net_QuasiStatic::net_QuasiStatic(int number):nNet(number){}

net_QuasiStatic::~net_QuasiStatic(){}


void net_QuasiStatic::start
(
	lexer *p, 
	fdm *a, 
	ghostcell *pgc
)
{
	double val_x, val_y, val_z, val_mag, vel_mag;

    // Reconstruct current net
    buildNet(p);
            
    // Calculating velocities at knots
    updateVel(p, a, pgc, 0);
    updateVel(p, a, pgc, 1);	
    updateVel(p, a, pgc, 2);


	// Solving the system of equations
	for (int i = 1; i < 500; i++)
	{	
		// Calculating force directions and coefficients
		for (int j = 0; j < nf; j++)
		{
            
			val_x = (v[Pi[j]][0] + v[Ni[j]][0])/2.0;
			val_y = (v[Pi[j]][1] + v[Ni[j]][1])/2.0;
			val_z = (v[Pi[j]][2] + v[Ni[j]][2])/2.0;
			
			vel_mag = sqrt(val_x*val_x + val_y*val_y + val_z*val_z);

			double vt = val_x*fi[j][0] + val_y*fi[j][1] + val_z*fi[j][2];
			double vt_mag = fabs(vt);

			e_d[j][0] = vt_mag*vt*fi[j][0];
			e_d[j][1] = vt_mag*vt*fi[j][1];
			e_d[j][2] = vt_mag*vt*fi[j][2];

			double vvtt_x = val_x - vt*fi[j][0];
			double vvtt_y = val_y - vt*fi[j][1];
			double vvtt_z = val_z - vt*fi[j][2];
			
			vt_mag = sqrt(vvtt_x*vvtt_x + vvtt_y*vvtt_y + vvtt_z*vvtt_z);
			
			e_l[j][0] = vt_mag*vvtt_x;
			e_l[j][1] = vt_mag*vvtt_y;
			e_l[j][2] = vt_mag*vvtt_z;


            /*
			e[j][0] = val_x/vel_mag;	
			e[j][1] = val_y/vel_mag;	
			e[j][2] = val_z/vel_mag;

			// Drag force
			e_d[j][0] = -e[j][0];
			e_d[j][1] = -e[j][1];
			e_d[j][2] = -e[j][2];			
			
			// Shear force
			val_x = fi[j][1]*e[j][2] - fi[j][2]*e[j][1];
			val_y = fi[j][2]*e[j][0] - fi[j][0]*e[j][2];
			val_z = fi[j][0]*e[j][1] - fi[j][1]*e[j][0];
				
			val_mag = sqrt(val_x*val_x + val_y*val_y + val_z*val_z);
				
			e_q[j][0] = val_x/val_mag;
			e_q[j][1] = val_y/val_mag;
			e_q[j][2] = val_z/val_mag;
	
			//	Lift force	
			val_x = e_q[j][1]*e[j][2] - e_q[j][2]*e[j][1];
			val_y = e_q[j][2]*e[j][0] - e_q[j][0]*e[j][2];
			val_z = e_q[j][0]*e[j][1] - e_q[j][1]*e[j][0];
			
			val_mag = sqrt(val_x*val_x + val_y*val_y + val_z*val_z);
				
			e_l[j][0] = val_x/val_mag;
			e_l[j][1] = val_y/val_mag;
			e_l[j][2] = val_z/val_mag;

			// Calculating force coefficients
			double arg = 
				(e[j][0]*fi[j][0] + e[j][1]*fi[j][1] + e[j][2]*fi[j][2])
				/(vel_mag*sqrt(fi[j][0]*fi[j][0] + fi[j][1]*fi[j][1] + fi[j][2]*fi[j][2]));
				
			if (arg > 1.0)
			{
				getC(acos(PI/2.0 - arg), c[j]);
			}
			else if (arg < 0.0)
			{
				if (fabs(arg) > 1.0)
				{
					getC(PI - acos(PI/2.0 - fabs(arg)), c[j]);
				}
				else
				{
					getC(PI - acos(fabs(arg)), c[j]);
				}
			}
			else
			{
				getC(acos(arg), c[j]);
			}
            */
		}		
		
		// Fill right-hand side
		int index = 0;
		bool bk;
		
		for (int j = 0; j < nK; j++)
		{
			bk = false;
		
			for (int k = 0; k < 2*n+2*m; k++)
			{
				if (j == Pb[k] || j == Nb[k])
				{
					bk = true;
					break;
				}
			}
			
			if (bk == false)
			{
				Bh[index][0] = B[index][0];
				Bh[index][1] = B[index][1];
				Bh[index][2] = B[index][2];		

				for (int k = 0; k < 4; k++)
				{
					int nfKik = nfK[index][k]; 	

					double cdt = 1.0;
					double cdn = 0.1;
					
					Bh[index][0] += p->W1/2.0*d_c*l[nfKik]/2.0*(cdt*e_d[nfKik][0] + cdn*e_l[nfKik][0]);
					Bh[index][1] += p->W1/2.0*d_c*l[nfKik]/2.0*(cdt*e_d[nfKik][1] + cdn*e_l[nfKik][1]);
					Bh[index][2] += p->W1/2.0*d_c*l[nfKik]/2.0*(cdt*e_d[nfKik][2] + cdn*e_l[nfKik][2]);

/*
					Bh[index][0] += 0.5*p->W1*v[index][0]*v[index][0]*d_c*l[nfKik]/2.0*
						(
							 c[nfKik][0]*e_d[nfKik][0]
							+ c[nfKik][1]*e_q[nfKik][0]
							+ c[nfKik][2]*e_l[nfKik][0]
						);
					
					Bh[index][1] += 0.5*p->W1*v[index][1]*v[index][1]*d_c*l[nfKik]/2.0*
						(
							 c[nfKik][0]*e_d[nfKik][1]
							+ c[nfKik][1]*e_q[nfKik][1]
							+ c[nfKik][2]*e_l[nfKik][1]
						);

					Bh[index][2] += 0.5*p->W1*v[index][2]*v[index][2]*d_c*l[nfKik]/2.0*
						(
							 c[nfKik][0]*e_d[nfKik][2]
							+ c[nfKik][1]*e_q[nfKik][2]
							+ c[nfKik][2]*e_l[nfKik][2]
						);
*/
				}
							
				index++;
			}
		} 

/*
MatrixXd A_eigen = MatrixXd::Zero(nf,nf);
MatrixXd B_eigen = MatrixXd::Zero(nf,3);
MatrixXd x_eigen = MatrixXd::Zero(nf,3);

double starttime1=pgc->timer();
for (int aa=0; aa<nf; aa++)
{
	B_eigen(aa,0) = Bh[aa][0];
	B_eigen(aa,1) = Bh[aa][1];
	B_eigen(aa,2) = Bh[aa][2];
	
	for (int bb=0; bb<nf; bb++)
	{
		A_eigen(aa,bb) = A[aa][bb];
	}
}
x_eigen = A_eigen.lu().solve(B_eigen);
*/

		// Solving the system
		solveGauss(p, A, Bh, fi);

		// Check error norm
		double norm, maxnorm;
			
		maxnorm = 0.0;
		for (int j = 0; j < nf; j++)
		{
			norm = 
				sqrt
				(
					 fi[j][0]*fi[j][0] + fi[j][1]*fi[j][1] + fi[j][2]*fi[j][2]
				);
			
			if (norm > maxnorm) maxnorm = norm;
		}
			
		if (fabs(maxnorm - 1.0) < 1e-3)
		{
			if (p->mpirank == 0)
			{
				cout<<"Number of iterations = "<<i<<setprecision(5)<<" with error = "<<fabs(maxnorm - 1.0)<<endl;
			}
			break;
		}
			
		// Otherwise correct system
		for (int j=0; j<nf; j++)	
		{
			norm = 
				sqrt
				(
					 fi[j][0]*fi[j][0] 
					+ fi[j][1]*fi[j][1] 
					+ fi[j][2]*fi[j][2]
				);
				
			fi[j][0] /= norm;
			fi[j][1] /= norm;
			fi[j][2] /= norm;
        
			for (int k = 0; k < niK; k++) 
			{
				A[k][j] *= norm;
			}
		}
		
		// Correct length of bars
		updateLength();
		
		// Reconstruct current net
		buildNet(p);
	}
	

	// Calculate horizontal and vertical reaction forces at mooring point
	//Xme_ = -fabs(A[sigma-1][sigma])*f[sigma][0];
	//Yme_ = -fabs(A[sigma-1][sigma])*f[sigma][1];
	//Zme_ = -fabs(A[sigma-1][sigma])*f[sigma][2];
		
	// Plotting net	
	print(p);			
}


void net_QuasiStatic::getC(double theta, double*& c_)
{
	theta *= 180/PI;
	
	c_[0] = max(0.0, -0.000000016676062*pow(theta,4.0) + 0.000001614232301*pow(theta,3.0) + 0.000132218750548*pow(theta,2.0) - 0.001335313221272*theta + 0.023082880162714);
	c_[1] = max(0.0, 0.000000019889794877*pow(theta,4.0) - 0.000005200545357615*pow(theta,3.0) + 0.000309491763289986*pow(theta,2.0) + 0.000866335727454720*theta - 0.000960622162325914);
	c_[2] = max(0.0, -0.000000028334918*pow(theta,4.0) + 0.000000348243732*pow(theta,3.0) + 0.000255407378882*pow(theta,2.0) - 0.006366497021954*theta + 0.049082170481768);
}


void net_QuasiStatic::updateLength()
{
	double T;
	
	for (int j = 0; j < nf; j++)	
	{
		T = 0.0;
		
		for (int i = 0; i < niK; i++)	
		{	
			T = MAX(fabs(A[i][j]),T); 
		}
		
		l[j] = l0[j]*(1.0 + T/EA);
		
		for (int i = niK; i < nf; i++)	
		{	
			if (A[i][j] > 0.0)
			{
				A[i][j] = l[j];
			}
			else if (A[i][j] < 0.0)
			{
				A[i][j] = -l[j];
			}
		}
	}
}


void net_QuasiStatic::updateVel(lexer *p, fdm *a, ghostcell *pgc, int cmp)
{
	int *recVel, *count;
	
	p->Iarray(count,p->mpi_size);
	p->Iarray(recVel, nK);
	
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
			if (cmp == 0)
			{
				v[i][cmp] = p->ccipol1_a(a->u,K_[i][0],K_[i][1],K_[i][2]);
			}
			else if (cmp == 1)
			{
				v[i][cmp] = p->ccipol2_a(a->v,K_[i][0],K_[i][1],K_[i][2]);
			}
			else
			{
				v[i][cmp] = p->ccipol3_a(a->w,K_[i][0],K_[i][1],K_[i][2]);
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
					K_[i][0] >= xstart[j] && K_[i][0] < xend[j] &&
					K_[i][1] >= ystart[j] && K_[i][1] < yend[j] &&
					K_[i][2] >= zstart[j] && K_[i][2] < zend[j]
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
	for (int i = 0; i < nK; i++)
	{
		if (recVel[i] == -1)
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
		
	for (int i = 0; i < nK; i++)
	{
		for (int j = 0; j < p->mpi_size; j++)
		{			
			if (recVel[i] == j)
			{		
				v[i][cmp] = recvVel[j][count[j]];
				count[j]++;
			}
		}
	}
	
	for (int i = 0; i < nK; i++)
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
	p->del_Iarray(recVel, nK);
}


void net_QuasiStatic::solveGauss
(
	lexer *p,
	double**& A, 
	double**& B,
    double**& f
) 
{
// n -> nf
	for (int i = 0; i < nf; i++)
	{
		for (int j = 0; j < nf; j++)
		{
			A_[i][j] = A[i][j];
		}
		
		B_[i][0] = B[i][0];
		B_[i][1] = B[i][1];
		B_[i][2] = B[i][2];
	}

    for (int i=0; i<nf; i++) 
	{
        // Search for maximum in this column
        double maxEl = fabs(A_[i][i]);
        int maxRow = i;
		
        for (int k=i+1; k<nf; k++) 
		{
            if (fabs(A_[k][i]) > maxEl) 
			{
                maxEl = fabs(A_[k][i]);
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        for (int k=i; k<nf; k++)
		{
            double tmp = A_[maxRow][k];
			
            A_[maxRow][k] = A_[i][k];
            A_[i][k] = tmp;
        }
        for (int k=0; k<3; k++)
		{
            double tmp = B_[maxRow][k];
			
            B_[maxRow][k] = B_[i][k];
            B_[i][k] = tmp;
        }
		
		
        // Make all rows below this one 0 in current column
        for (int k=i+1; k<nf; k++) 
		{
            double c = -A_[k][i]/A_[i][i];
			
            for (int j=i; j<nf; j++)
			{
                if (i==j) 
				{
                    A_[k][j] = 0;
				} 
				else 
				{
                    A_[k][j] += c * A_[i][j];
                }
            }
            for (int j=0; j<3; j++)
			{
				B_[k][j] += c * B_[i][j];
            }			
        }
    }

    // Solve equation A_*f=B_ for an upper triangular matrix A_
	
    for (int i=nf-1; i>=0; i--) 
	{
		for (int j=0; j<3; j++)
		{
			f[i][j] = B_[i][j]/A_[i][i];
		}			
		
        for (int k=i-1; k>=0; k--) 
		{
			for (int j=0; j<3; j++)
			{
				B_[k][j] -= A_[k][i] * f[i][j];
			}			
        }
    }
}


void net_QuasiStatic::netForces
(
	double& Xme, double& Yme, double& Zme
)
{
	Xme = 0.0;
	Yme = 0.0;
	Zme = 0.0;
}
