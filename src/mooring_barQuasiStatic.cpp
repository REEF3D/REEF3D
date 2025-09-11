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
	
	// Calculate horizontal and vertical reaction forces at mooring point
	Xme_ = -fabs(A[sigma-1][sigma])*f[sigma][0];
	Yme_ = -fabs(A[sigma-1][sigma])*f[sigma][1];
	Zme_ = -fabs(A[sigma-1][sigma])*f[sigma][2];
		
	// Plotting mooring line	
	print(p,pgc);	
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

void mooring_barQuasiStatic::mooringForces(double& Xme, double& Yme, double& Zme)
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
