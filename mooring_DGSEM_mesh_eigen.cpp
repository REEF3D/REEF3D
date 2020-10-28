/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2020 Tobias Martin

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
#include <Eigen/Dense>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd;

void mooring_DGSEM::facenodegen(lexer *p, double xS, double xE)
{
	p->Darray(vx, H+1);

	for (int i = 0; i < (H+1); i++)
	{
		vx[i] = xS + i*(xE-xS)/H;
	}
}

void mooring_DGSEM::getqp(lexer *p)
{
	p->Darray(r, P+1);

	r[0] = -1.0;
	r[P] = 1.0;
		
	if (P == 2)
	{
		r[0] = -1.0;
		r[1] = 0.0;
		r[2] = 1.0;
	}
	else if (P > 2)
	{
		VectorXd h1(P-1);
		for (int i = 0; i < h1.size(); i++)
		{
			h1(i) = 2.0 + 2.0*i;
		}
		
		MatrixXd J = MatrixXd::Zero(P-1,P-1);	
		for (int i = 0; i < (P-2); i++)
		{		
			J(i, i+1) = 
				2.0/(h1(i) + 2.0)
				*sqrt((i+1)*(i+3)*(i+2)*(i+2)/((h1(i) + 1.0)*(h1(i) + 3.0)));
		}
		
		MatrixXd Jtranspose = J.transpose();
		J = J + Jtranspose;
			
		// Calculate eigenvalues of J and sort them
		Eigen::EigenSolver<MatrixXd> eigensolver(J);
		for (int i = 0; i < (P-1); i++)
		{
			r[i+1] = eigensolver.eigenvalues().col(0)[i].real();
		}

		double temp = 0.0;
		for(int i = 0; i < (P+1); i++)
		{		
			for(int j = (i+1); j < (P+1); j++)
			{
				if(r[i] > r[j])
				{
					temp = r[i];
					r[i] = r[j];
					r[j] = temp;
				}
			}
		}
		
		J.resize(0,0);
		h1.resize(0);
	}
}


void mooring_DGSEM::getDr(lexer *p)
{
	p->Darray(Dr,P+1,P+1);
	p->Darray(V,P+1,P+1);
	p->Darray(invV,P+1,P+1);

	MatrixXd V_ = MatrixXd::Zero(P+1, P+1);
	MatrixXd Vr_ = MatrixXd::Zero(P+1, P+1);
	MatrixXd Dr_ = MatrixXd::Zero(P+1, P+1);
		
	double *h;
	p->Darray(h, P+1);
		
	// Calculate V
	for (int i = 1; i < (P + 2); i++)
	{			
		jacobiP(0, 0, i-1, h);
			
		for (int j = 0; j < (P + 1); j++)
		{
			V_(j,i-1) = h[j]; 
		}
	}
		
	// Calculate Vr
	for (int i = 0; i < (P + 1); i++)
	{
		if (i == 0)
		{
			for (int j = 0; j < (P + 1); j++)
			{
				Vr_(j,i) = 0.0; 
			}
		}
		else
		{
			jacobiP(1.0, 1.0, i-1, h);

			for (int j = 0; j < (P + 1); j++)
			{
				Vr_(j,i) = sqrt(i*(i+1))*h[j]; 
			}	
		}
	}
			
	// Dr = Vr/V;
	Dr_ = Vr_*V_.inverse();
	
	
	// Communicate to solver
	MatrixXd Vinv = V_.inverse();
	
	for (int i = 0; i < (P + 1); i++)
	{
		for (int j = 0; j < (P + 1); j++)
		{
			V[i][j] = V_(i,j);
			invV[i][j] = Vinv(i,j);
		}
	}

	for (int i = 0; i < (P + 1); i++)
	{
		for (int j = 0; j < (P + 1); j++)
		{
			Dr[i][j] = Dr_(i,j);
		}
	}
		
	// Delete matrices
	V_.resize(0,0);
	Vr_.resize(0,0);
	Dr_.resize(0,0);
	p->del_Darray(h, P+1);
}


void mooring_DGSEM::jacobiP(double alpha, double beta, int N, double *& h)
{
	MatrixXd Pl = MatrixXd::Zero(N+1, P+1);

	double gamma0 = 
		pow(2.0,alpha+beta+1.0)/(alpha+beta+1.0)*tgamma(alpha+1.0)
		*tgamma(beta+1.0)/tgamma(alpha+beta+1.0);
		
	for (int i = 0; i < (P + 1); i++)
	{
		Pl(0,i) = 1.0/sqrt(gamma0);
	}
			
	if (N == 0) 
	{		
		for (int i = 0; i < (P + 1); i++)
		{
			h[i] = Pl(0,i);
		}
	}
	else
	{
		double gamma1 = (alpha+1.0)*(beta+1.0)/(alpha+beta+3.0)*gamma0;
				
		for (int i = 0; i < (P + 1); i++)
		{
			Pl(1,i) = ((alpha+beta+2.0)*r[i]/2.0 + (alpha-beta)/2.0)/sqrt(gamma1);
		}							

		if (N == 1) 
		{
			for (int i = 0; i < (P + 1); i++)
			{
				h[i] = Pl(1,i);
			}
		}	
		else
		{
			double aold = 2.0/(2.0+alpha+beta)*sqrt((alpha+1.0)*(beta+1.0)/(alpha+beta+3.0));
					
			for (int i = 1; i < (N); i++)
			{	
				double h1 = 2*i + alpha + beta;
						
				double anew = 2/(h1+2)*sqrt( (i+1)*(i+1+alpha+beta)*(i+1+alpha)*(i+1+beta)/(h1+1)/(h1+3));
				double bnew = - (alpha*alpha-beta*beta)/h1/(h1+2);
	
				for (int j = 0; j < (P+1); j++)
				{					
					Pl(i+1,j) = 1.0/anew*(-aold*Pl(i-1,j) + r[j]*Pl(i,j));
				}
  
				aold =anew;
			}

			for (int i = 0; i < (P + 1); i++)
			{
				h[i] = Pl(N,i);
			}					
		}
	}
			
	Pl.resize(0,0);	
}


void mooring_DGSEM::getsurfInt(lexer *p)
{
	p->Darray(sInt, P+1, 2);

	MatrixXd F_ = MatrixXd::Zero(P+1, 2);
	MatrixXd V_ = MatrixXd::Zero(P+1, P+1);
	MatrixXd sInt_ = MatrixXd::Zero(P+1, 2);

	for (int i = 0; i < (P + 1); i++)
	{
		for (int j = 0; j < (P + 1); j++)
		{
			V_(i,j) = V[i][j];
		}
	}
	F_(0,0) = 1.0;
	F_(P,1) = 1.0;
	
	sInt_ = V_*(V_.transpose()*F_);
	
	for (int i = 0; i < (P + 1); i++)
	{
		sInt[i][0] = sInt_(i,0);	
		sInt[i][1] = sInt_(i,1);	
	}

	F_.resize(0,0);
	V_.resize(0,0);
	sInt_.resize(0,0);
}


void mooring_DGSEM::getV(lexer *p){}


void mooring_DGSEM::nodegen(lexer *p)
{
	p->Darray(x,H,P+1);

	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < (P+1); j++)
		{
			x[i][j] = vx[i] + 0.5*(r[j] + 1.0)*(vx[i+1] - vx[i]);
		}
	}

	// Scale factor assuming constant delta H
	rx = 0.0;
	for (int j = 0; j < (P+1); j++)
	{	
		rx += Dr[0][j]*x[0][j];
	}
	rx = 1.0/(rx);
}
