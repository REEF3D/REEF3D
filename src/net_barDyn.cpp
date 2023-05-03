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

#include"net_barDyn.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"reinidisc_fsf.h"	

net_barDyn::net_barDyn(int number, lexer *p):nNet(number){}

net_barDyn::~net_barDyn(){}


void net_barDyn::initialize(lexer *p, fdm *a, ghostcell *pgc)
{    
    //- Initialise net model
    if (p->X320_type[nNet] == 12)   
    {
        cyl_ini(p,a,pgc);
 
        buildNet_cyl(p); 
    }
    else if (p->X320_type[nNet] == 13)   
    {
        wall_ini(p,a,pgc);
        
        buildNet_wall(p);    
    } 
    else if (p->X320_type[nNet] == 14)   
    {
        cone_ini(p,a,pgc);
        
        buildNet_wall(p);    
    } 

    //- Initialise old variables   
    xnn_ = x_;    
    xn_ = x_;
    xdotnn_ = xdot_;    
    xdotn_ = xdot_;
    
    if (p->X325_dt == 0.0)
    {   
        dtnn_ = p->dt;
        dtn_ = p->dt;
        dt_ = p->dt;
    }
    else
    {
        dtnn_ = p->X325_dt;
        dtn_ = p->X325_dt;
        dt_ = p->X325_dt;
    }
    t_net_n = 0.0;
    t_net = 0.0;

    //- Initialise printing
    printtime = 0.0;
    print(p);
}


void net_barDyn::start
(
	lexer *p, 
	fdm *a, 
	ghostcell *pgc,
    double alpha,
    Eigen::Matrix3d quatRotMat
)
{
    double starttime1 = pgc->timer();    

	//- Set net time step
	double phi = 0.0;
	t_net_n = t_net;
	t_net = phi*p->simtime + (1.0 - phi)*(p->simtime + alpha*p->dt);
	double dtm = t_net - t_net_n;
	
    dt_ = p->X325_dt > 0.0 ? min(dtm, p->X325_dt) : dtm;

	//- Start loop
    int loops = ceil(dtm/dt_);
    if (dt_ == 0.0) loops = 0;
    dt_ = dtm/loops;

	Eigen::VectorXi convIt(loops);
    
    for (int loop = 0; loop < loops; loop++)
    {
        convIt(loop) = loop;

        startLoop(p,a,pgc,convIt(loop));
    }

    //- Coupling forces for vrans model
    vransCoupling(p,a,pgc);

	//- Build and save net
	print(p);	

    //- Print output
    double endtime1 = pgc->timer() - starttime1; 
    if (p->mpirank == 0 && convIt.maxCoeff() < 10)
    {
        cout<<"Net converged within "<<convIt.maxCoeff()<<" iterations max, "<<loops<<" steps and "<<endtime1<<" s"<<endl;
    }
    if (p->mpirank == 0 && convIt.maxCoeff() >= 10)
    {
        cout<<"Net diverged. Adjust X 325 accordingly!"<<endl; 
    }
}


void net_barDyn::startLoop
(
	lexer *p, 
	fdm *a, 
	ghostcell *pgc,
    int& iter
)
{
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

    //- Get time weights for finite differences
    coeffs_ = timeWeight(p);

    //- Calculate force vector
    getForces(p);

    //- Solve dynamics
    if (p->count < 1 && iter == 0)
    {
        //- Solve linear system
        
        // Fill system of equations in A_
        fillLinSystem(p, a, pgc);

        // Fill RHS in B_
        fillLinRhs(p, a, pgc);

        // Solve system for tension forces
        T_ = A_.partialPivLu().solve(B_.transpose());
	    limitTension();

        iter = 1;
    }
    else
    {
        //- Solve non-linear system
        
        double norm_error;
       
        T_backup = T_;
        iter = 0;

        for (int it = 0; it < 10; it++)
        {
            // Fill Jacobian and invert
            fillNonLinSystem(p, a, pgc); 

            Eigen::PartialPivLU<MatrixXd> inv(A_);
            
            // Fill non-linear function
            fillNonLinRhs(p, a, pgc);

	        // Store tension forces
	        T_old = T_;

            // Solve system for intermediate tension forces
            T_ -= inv.solve(B_.transpose());
	        limitTension();

	        // Accelerated Newton step
	        fillNonLinRhs(p, a, pgc);
	        
            T_ -= inv.solve(B_.transpose());
            limitTension();

            // Check convergence
            norm_error = (T_ - T_old).norm();

            iter++;

            if (norm_error < 1e-10) 
	        { 
		        break;
            }
        }

        if (iter >= 9)
        {
            T_ = T_backup;
        }
    }
    
    //- Calculate accelerations
    updateAcc(p, a, pgc); 


    //- Advance velocitie
    xdot_ = 
        1.0/coeffs_(0)*
        (
            xdotdot_ - coeffs_(1)*xdot_ - coeffs_(2)*xdotn_ - coeffs_(3)*xdotnn_
        );

    if (p->X320 == 13)  // 2D solution
    {
        xdot_.col(1) *= 0.0;
    }

    /*MatrixXd xdot_new(nK,3); 
    
    xdot_new = 
        1.0/coeffs_(0)*
        (
            xdotdot_ - coeffs_(1)*xdot_ - coeffs_(2)*xdotn_ - coeffs_(3)*xdotnn_
        );
    */

    // Under-relaxation for velocities below rigid top knots
/*    xdot_new.block(nK-niK,0,niK,1) = xdot_.block(nK-niK,0,niK,1) + p->X325_relX*(xdot_new.block(nK-niK,0,niK,1) - xdot_.block(nK-niK,0,niK,1));
    xdot_new.block(nK-niK,1,niK,1) = xdot_.block(nK-niK,1,niK,1) + p->X325_relY*(xdot_new.block(nK-niK,1,niK,1) - xdot_.block(nK-niK,1,niK,1));
    xdot_new.block(nK-niK,2,niK,1) = xdot_.block(nK-niK,2,niK,1) + p->X325_relZ*(xdot_new.block(nK-niK,2,niK,1) - xdot_.block(nK-niK,2,niK,1));
*/

    //- Advance position
    x_ =         
        1.0/coeffs_(0)*
        (
            xdot_ - coeffs_(1)*x_ - coeffs_(2)*xn_ - coeffs_(3)*xnn_
        );
    /*
    x_new =         
        1.0/coeffs_(0)*
        (
            xdot_new - coeffs_(1)*x_ - coeffs_(2)*xn_ - coeffs_(3)*xnn_
        );
    */

    //- Save old velocity and position vectors
    xnn_ = xn_;    
    xn_ = x_;
    xdotnn_ = xdotn_;    
    xdotn_ = xdot_;
    //xdot_ = xdot_new;
    
    //- Save old time steps
    dtnn_ = dtn_;
    dtn_ = dt_;
    
}


void net_barDyn::limitTension()
{ 
    // Avoid unphysical tension forces
    for (int barI = 0; barI < nK; barI++)
    {
        T_(barI) = max(T_(barI), 0.0);
    }
}


void net_barDyn::updateAcc(lexer *p, fdm *a, ghostcell *pgc)
{
    Vector3d T_knot, x_ij;
    int knotJ, barI;
    double l_ij;

    // Move top as rigid body
    updateTopAcc(p);

    // Calculate acceleration from tension forces
    for (int knotI = 0; knotI < nK; knotI++)
    {
        T_knot << 0.0, 0.0, 0.0;
        
        if (knotI >= nfK[0][0]) // then inner knot
        {
            xdotdot_.row(knotI) *= 0.0; 
            
            for (int k = 1; k < 5; k++)
            {
                barI = nfK[knotI - nfK[0][0]][k];

                if (barI != -1)
                {
                    // Find bar vector
                    if (Pi[barI] == knotI)
                    {
                        knotJ = Ni[barI];
                    }
                    else
                    {
                        knotJ = Pi[barI];
                    }

                    x_ij = x_.row(knotJ) - x_.row(knotI);

                    l_ij = x_ij.norm();

                    T_knot += T_(barI)*x_ij/l_ij; 
                }
            }

            xdotdot_.row(knotI) = 
                (forces_knot.row(knotI) + T_knot)
                /
                (mass_knot(knotI) + added_mass(knotI));
        }
        else    // rigid motion of top knots
        {
            xdotdot_.row(knotI) = top_xdotdot_.row(knotI);
        }
    } 
}


void net_barDyn::updateTopAcc(lexer *p)
{
    // Update top knot velocities
    for (int i = 0; i < nbK; i++)
    {
        top_xdot_(i,0) = p->ufbi + (x_(i,2) - p->zg)*p->qfbi - (x_(i,1) - p->yg)*p->rfbi;
        top_xdot_(i,1) = p->vfbi + (x_(i,0) - p->xg)*p->rfbi - (x_(i,2) - p->zg)*p->pfbi;
        top_xdot_(i,2) = p->wfbi + (x_(i,1) - p->yg)*p->pfbi - (x_(i,0) - p->xg)*p->qfbi;
    }

    // Calculate top knot acceleration
    top_xdotdot_ = 
          coeffs_(0)*top_xdot_ + coeffs_(1)*xdot_.block(0,0,nbK,3)
        + coeffs_(2)*xdotn_.block(0,0,nbK,3) + coeffs_(3)*xdotnn_.block(0,0,nbK,3);
}


void net_barDyn::updateField(lexer *p, fdm *a, ghostcell *pgc, int cmp)
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


Eigen::VectorXd net_barDyn::timeWeight(lexer* p)
{	
    // 3rd-order finite difference weights for first derivative and varying time step
    
	double c2, c3, c5;
	int mn;	

	int nd = 4;    
    double c1 = 1.0;
    double c4 = 0.0;

    if (p->count == 1)
    {
        dtn_ = dt_;
        dtnn_ = dt_;
    }
    
    VectorXd ti(nd);

	ti(0) = 0.0;
    ti(1) = ti(0) - dt_;
    ti(2) = ti(1) - dtn_;
    ti(3) = ti(2) - dtnn_;
    
    MatrixXd coeff = MatrixXd::Zero(nd,nd);    
            
    coeff(0,0) = 1.0;

    for(int r = 1; r < nd; ++r)
    {
        int mn = MIN(r,1);
        
        c2 = 1.0;
        c5 = c4;
        c4 = ti(r) - ti(0);

        for(int s = 0; s < r; ++s)
        {
            c3 = ti(r) - ti(s);
            c2 *= c3;
                
            if(s == r-1) 
            {
                for(int t = mn; t >= 1; t--)
                {
                    coeff(r,t) = c1*(double(t)*coeff(r-1,t-1) - c5*coeff(r-1,t))/c2;
                }
         
                coeff(r,0) = -c1*c5*coeff(r-1,0)/c2;
            }
                
            for(int t = mn; t >= 1; t--)
            {
                coeff(s,t) = (c4*coeff(s,t) - double(t)*coeff(s,t-1))/c3;
            }
           
            coeff(s,0) *= c4/c3;
        }
        c1 = c2;
    }
    
    return coeff.col(1);    
}
