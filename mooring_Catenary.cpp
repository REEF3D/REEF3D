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

#include"mooring_Catenary.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

mooring_Catenary::mooring_Catenary(int number):line(number){}

mooring_Catenary::~mooring_Catenary(){}

void mooring_Catenary::start(lexer *p, fdm *a, ghostcell *pgc)
{
	// Calculate distances between start and mooring points
	dx = p->X311_xe[line] - p->X311_xs[line];			
	dy = p->X311_ye[line] - p->X311_ys[line];				
	dz = p->X311_ze[line] - p->X311_zs[line];	

	dxy = sqrt(dx*dx+dy*dy);			


	// Calculating forces at mooring point
	if (EA == 0.0)
	{
		Fh_0= w + 1.0;
		FH = w;
		while (fabs(FH - Fh_0) > 1e-3)
		{
			Fh_0 = FH;

			f_Fh = 
				L 
				- dz*sqrt(1.0 + 2.0*FH/(dz*w)) + FH/w*acosh(1.0 + w*dz/FH) 
				- dxy;
			ddf_Fh = 
				-dz/(sqrt(2.0 + (dz*w)/FH)*sqrt((dz*w)/FH)*FH) 
				- 1.0/(w*sqrt(1.0 + (2.0*FH)/(dz*w))) 
				+ (acosh(1.0 + (dz*w)/FH))/w;

			FH -= f_Fh/ddf_Fh;
		}
		
		Xme_ = FH*fabs(cos(atan(dy/dx)));
		Yme_ = FH*fabs(sin(atan(dy/dx)));

		lms = dz*sqrt(1.0 + 2.0*FH/(dz*w));
		Zme_ = w*lms;
	}
	else
	{
		double dxy_ = dxy;

		for (int loop = 0; loop < 1000; loop++)
		{
			double f1, f2, df1H, df1V, df2H, df2V;
			FH = 0.001;
			FV = 0.001;
			double FH_0 = FH + 1.0;
			double FV_0 = FV + 1.0;	
			
			while (fabs(FH - FH_0) > 1e-5 && fabs(FV - FV_0) > 1e-5)
			{
				FH_0 = FH;
				FV_0 = FV;
		
				f1 = L - FV/w + FH/w*log(FV/FH + sqrt(1+(FV/FH)*(FV/FH))) - dxy_;
				f2 = FH/w*(sqrt(1+(FV/FH)*(FV/FH))-1) + FV*FV/(2*EA*w) - dz;
			   
				df1H = 1/w*log(FV/FH + sqrt(1+(FV/FH)*(FV/FH))) + (FH*(-FV/(FH*FH) - FV*FV/(FH*FH*FH*sqrt(FV*FV/(FH*FH)+1))))/(w*(FV/FH + sqrt(1+(FV/FH)*(FV/FH))));
				df1V = -1/w + 1/(w*sqrt((FV/FH)*(FV/FH)+1));
				df2H = 1/w*(sqrt(1+(FV/FH)*(FV/FH))-1) - FV*FV/(FH*FH*w*sqrt((FV/FH)*(FV/FH)+1));
				df2V = FV/(EA*w) + FV/(FH*w*sqrt((FV/FH)*(FV/FH)+1));

				A[0][0] = df1H;
				A[0][1] = df1V;
				A[1][0] = df2H;
				A[1][1] = df2V;
				B[0] = -f1;
				B[1] = -f2;
		
				for (int i=0; i<2; i++) 
				{
					// Search for maximum in this column
					double maxEl = fabs(A[i][i]);
					int maxRow = i;
					
					for (int k=i+1; k<2; k++) 
					{
						if (fabs(A[k][i]) > maxEl) 
						{
							maxEl = fabs(A[k][i]);
							maxRow = k;
						}
					}

					// Swap maximum row with current row (column by column)
					for (int k=i; k<2; k++)
					{
						double tmp = A[maxRow][k];
						
						A[maxRow][k] = A[i][k];
						A[i][k] = tmp;
					}
					double tmp = B[maxRow];
						
					B[maxRow] = B[i];
					B[i] = tmp;
					
					// Make all rows below this one 0 in current column
					for (int k=i+1; k<2; k++) 
					{
						double c = -A[k][i]/A[i][i];
						
						for (int j=i; j<2; j++)
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

						B[k] += c * B[i];
					}
				}

				// Solve equation Ax=b for an upper triangular matrix A
				for (int i=1; i>=0; i--) 
				{
					F[i] = B[i]/A[i][i];
						
					for (int k=i-1; k>=0; k--) 
					{
						B[k] -= A[k][i] * F[i];		
					}
				}

				FH = FH + F[0];
				FV = FV + F[1];
			}
		
			Xme_ = FH*fabs(cos(atan(dy/dx)));
			Yme_ = FH*fabs(sin(atan(dy/dx)));
			Zme_ = FV;	
			lms = L - FV/w;
			
			// Check convergence
			buildLine(p);
			
			if (fabs(dxy - x[H-1]) > 0.01)
			{
				if (x[H-1] - dxy > 0.01)
				{
					dxy_ -= 0.01;
				}
				else
				{
					dxy_ += 0.01;
				}
			}
			else
			{
				break;
			}
		}
	}
	

	// Reaction forces at mooring points	
	if (dx > 0)	
	{
		Xme_ *= -1.0;
	}
	if (dy > 0)	
	{
		Yme_ *= -1.0;
	}
	if (dz > 0)	
	{
		Zme_ *= -1.0;
	}
		
	// Print mooring line
	print(p);
}


void mooring_Catenary::mooringForces
(
	double& Xme, double& Yme, double& Zme
)
{
	Xme = Xme_; 
	Yme = Yme_;
	Zme = Zme_;
}


void mooring_Catenary::iniDyn(lexer *p, fdm *a, ghostcell *pgc, double& FH_, double& FV_)
{
	double rho_f = 1000.0;
	
	rho_c = p->X311_rho_c[line];
	w = p->X311_w[line]*9.81*(rho_c - rho_f)/rho_c;
	L = p->X311_l[line];
	H = p->X311_H[line];
	EA = p->X311_EA[line];
	
	xs = p->X311_xs[line];
	ys = p->X311_ys[line];
	zs = p->X311_zs[line];
	
	p->Darray(x,H); 
	p->Darray(y,H);
	p->Darray(z,H); 
	p->Darray(T,H);
	p->Darray(B,2);
	p->Darray(F,2);
	p->Darray(A,2,2);

	printtime = 0.0;

	start(p, a, pgc);
	
	FH_ = FH;
	FV_ = FV;
}
