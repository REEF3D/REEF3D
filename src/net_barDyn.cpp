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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"net_barDyn.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm_nhf.h"
#include"ghostcell.h"	

net_barDyn::net_barDyn(int number, lexer *p):nNet(number)
{
}

net_barDyn::~net_barDyn()
{
}

void net_barDyn::start_cfd(lexer *p, fdm *a, ghostcell *pgc, double alpha, Eigen::Matrix3d quatRotMat)
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
    if (dt_==0.0) loops = 0;
    dt_ = dtm/loops;

	Eigen::VectorXi convIt(loops);
    
    for (int loop = 0; loop < loops; loop++)
    {
        convIt(loop) = loop;
        
        update_velocity_cfd(p,a,pgc);
        startLoop(p,pgc,convIt(loop));
    }

    //- Coupling forces for vrans model
    coupling_dlm_cfd(p,a,pgc);

	//- Build and save net
	print(p);	
}

void net_barDyn::start_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha, Eigen::Matrix3d quatRotMat)
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
    if (dt_==0.0) loops = 0;
    dt_ = dtm/loops;

	Eigen::VectorXi convIt(loops);
    
    for (int loop = 0; loop < loops; loop++)
    {
        convIt(loop) = loop;
        
        update_velocity_nhflow(p,d,pgc);
        startLoop(p,pgc,convIt(loop));
    }

    //- Coupling forces for vrans model
    coupling_dlm_nhflow(p,d,pgc);

	//- Build and save net
	print(p);	
}


void net_barDyn::startLoop(lexer *p, ghostcell *pgc, int& iter)
{
          
    //- Get time weights for finite differences
    coeffs_ = timeWeight(p);

    //- Calculate force vector
    getForces(p);

    //- Solve dynamics
    if (p->count < 1 && iter==0)
    {
        //- Solve linear system
        
        // Fill system of equations in A_
        fillLinSystem(p, pgc);

        // Fill RHS in B_
        fillLinRhs(p, pgc);

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
            fillNonLinSystem(p, pgc); 

            Eigen::PartialPivLU<MatrixXd> inv(A_);
            
            // Fill non-linear function
            fillNonLinRhs(p, pgc);

	        // Store tension forces
	        T_old = T_;

            // Solve system for intermediate tension forces
            T_ -= inv.solve(B_.transpose());
	        limitTension();

	        // Accelerated Newton step
	        fillNonLinRhs(p, pgc);
	        
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
    updateAcc(p, pgc); 


    //- Advance velocitie
    xdot_ = 
        1.0/coeffs_(0)*
        (
            xdotdot_ - coeffs_(1)*xdot_ - coeffs_(2)*xdotn_ - coeffs_(3)*xdotnn_
        );

    if (p->X320==13)  // 2D solution
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

Eigen::VectorXd net_barDyn::timeWeight(lexer* p)
{	
    // 3rd-order finite difference weights for first derivative and varying time step
    
	double c2, c3, c5;
	int mn;	

	int nd = 4;    
    double c1 = 1.0;
    double c4 = 0.0;

    if (p->count==1)
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
                
            if(s==r-1) 
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
