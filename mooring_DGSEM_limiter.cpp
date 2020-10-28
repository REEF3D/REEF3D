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

void mooring_DGSEM::limitFields
(
	lexer *p, 
	double **& r_x_, double **& r_y_, double **& r_z_, 
	double **& q_x_, double **& q_y_, double **& q_z_, 
	double **& v_x_, double **& v_y_, double **& v_z_
)
{
/*	startLimit(p, r_x_); 
	startLimit(p, r_y_); 
	startLimit(p, r_z_);
	
	startLimit(p, q_x_); 
	startLimit(p, q_y_); 
	startLimit(p, q_z_); 
	
	startLimit(p, v_x_); 
	startLimit(p, v_y_); 
	startLimit(p, v_z_);*/
}

	
void mooring_DGSEM::startLimit
(
	lexer *p, 
	double **& q
)
{
	double avg1, avg2;
	double **qlimit;
	double *qavg, *qlin, *lin_coeff;
	
	p->Darray(qlimit, H, P+1);
	p->Darray(qlin, P+1);
	p->Darray(lin_coeff, P+1);
	p->Darray(qavg, H+2);

	// Calculate average values in each cell plus ghostcells for +-1 operations
    for (int i = 0; i < H; i++)
    {
		for (int j = 0; j < (P+1); j++)
		{
			qlimit[i][j] = q[i][j];
			qavg[i+1] += invV[0][j]*q[i][j];
		}
		
		qavg[i+1] *= V[0][0];
    }
	qavg[0] = - qavg[1];														// set it rather as 2*a - qavg[1]
	qavg[H+1] = qavg[H];														// set it rather as 2*b - qavg[H]
	
		
	// Limiting process
    for (int i = 0; i < H; i++)
    {
		int iavg = i+1;
		
		avg1 = 
			qavg[iavg] - minMod
			(
				qavg[iavg] - q[i][0], 
				qavg[iavg] - qavg[iavg-1], 
				qavg[iavg+1] - qavg[iavg]
			);
			
		avg2 = 
			qavg[iavg] + minMod
			(
				q[i][P] - qavg[iavg], 
				qavg[iavg] - qavg[iavg-1], 
				qavg[iavg+1] - qavg[iavg]
			);
	
		// Cell needs limiting
		if((fabs(avg1 - q[i][0]) > 1.0e-8 || fabs(avg2 - q[i][P]) > 1.0e-8))
		{
			
		
			// Get linear coefficients in cell i
			for (int k = 0; k < (P+1); k++)
			{
				if (k < 2)
				{
					lin_coeff[k] = 0.0;
					
					for (int j = 0; j < (P+1); j++)
					{
						lin_coeff[k] += invV[k][j]*q[i][j];
					}
				}
				else
				{
					lin_coeff[k] = 0.0;
				}
			}
					
			// Reconstruct linear function qlin in cell i
			for (int k = 0; k < (P+1); k++)
			{
				qlin[k] = 0.0;
				
				for (int j = 0; j < (P+1); j++)
				{
					qlin[k] += V[k][j]*lin_coeff[j];
				}
			}			

			slopeLimit(p, i, qlimit[i], qlin, qavg);
		}
	}
	
	// Assign limited values
    for (int i = 0; i < H; i++)
    {
		for (int j = 0; j < (P+1); j++)
		{
			q[i][j] = qlimit[i][j];
			
		}
	}
	
	// Delete arrays
	p->del_Darray(qlimit, H, P+1);
	p->del_Darray(qlin, P+1);
	p->del_Darray(lin_coeff, P+1);
	p->del_Darray(qavg, H+2);
}


void mooring_DGSEM::slopeLimit
(
	lexer *p,
	int i,
	double*& qlimitI, 
	double*& qlin, 
	double*& qavg
)
{
	int iavg = i+1;
	double qdx = 0.0;
	double h = (x[i][P] - x[i][0]);
	
	// Calculate cell centre
	double x0 = x[i][0] + h/2.0;

	// Calculate limited coefficients
	for (int k = 0; k < (P+1); k++)
	{
		qdx = 0.0;
		
		for (int j = 0; j < (P+1); j++)
		{
			qdx += (Dr[k][j]*qlin[j]);
		}
		qdx *= 2.0/h;
		
		qlimitI[k] = 
			qavg[iavg] + (x[i][k] - x0)*minMod
			(
				qdx, 
				(qavg[iavg+1] - qavg[iavg])/h, 
				(qavg[iavg] - qavg[iavg-1])/h
			);			
	}
}


double mooring_DGSEM::minMod(const double& a, const double& b, const double& c)
{
	double s = 1.0/3.0*(SIGN(a) + SIGN(b) + SIGN(c));
	
	if (fabs(s) == 1)
	{
		s = s*MIN(MIN(fabs(a),fabs(b)),fabs(c));
	}
	else
	{
		s = 0.0;
	}
	
	return s;
}