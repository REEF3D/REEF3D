/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"VOF_PLIC.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"solver.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"fluid_update_vof.h"
#include"heat.h"
#include"hires.h"
#include"weno_hj.h"
#include"hric.h"

void VOF_PLIC::updateVOF(fdm* a, lexer* p, int sweep)
{
	if (sweep == 0)
	{
		LOOP
		{
			a->vof(i,j,k) = 
				max(0.0,min(vof3(i-1,j,k) + vof2(i,j,k) + vof1(i+1,j,k), 1.0));
		}
	}
	else if (sweep == 1)
	{
		LOOP
		{
			a->vof(i,j,k) = 
				max(0.0,min(vof3(i,j-1,k) + vof2(i,j,k) + vof1(i,j+1,k), 1.0));
		}
	}
	else
	{
		LOOP
		{			
			a->vof(i,j,k) = 
				max(0.0,min(vof3(i,j,k-1) + vof2(i,j,k) + vof1(i,j,k+1), 1.0));			
		}

	}
}


void VOF_PLIC::updateVolumeFraction
(
	fdm* a, 
	lexer* p, 
	const double Q1, 
	const double Q2, 
	int sweep
)
{
	//- Calculate volume entering, leaving and staying in the cell
	
	double m1 = max(Q1, 0.0);
	double m2 = 1.0 - m1 + min(0.0, Q2);

	if (sweep == 0)
	{
		if (Q1 < 0.0)
		{
			vof1(i, j, k) = 
				calcV
				(
					nx(i, j, k),
					ny(i, j, k),
					nz(i, j, k),
					alpha(i, j, k),
					Q1,
					-Q1
				);

		}
		if (Q2 > 0.0) 
		{
			vof3(i, j, k) = 
				calcV
				(
					nx(i, j, k),
					ny(i, j, k),
					nz(i, j, k),
					alpha(i, j, k),
					1.0,
					Q2
				);
		}
		
		vof2(i, j, k) = 
				calcV
				(
					nx(i, j, k),
					ny(i, j, k),
					nz(i, j, k),
					alpha(i, j, k),
					m1,
					m2
				);
	}
	else if (sweep == 1)
	{
		if (Q1 < 0.0)
		{
			vof1(i, j, k) = 
				calcV
				(
					ny(i, j, k),
					nz(i, j, k),
					nx(i, j, k),
					alpha(i, j, k),
					Q1,
					-Q1
				);
		}
		if (Q2 > 0.0) 
		{
			vof3(i, j, k) = 
				calcV
				(
					ny(i, j, k),
					nz(i, j, k),
					nx(i, j, k),
					alpha(i, j, k),
					1.0,
					Q2
				);
		}
		
		vof2(i, j, k) = 
			calcV
			(
				ny(i, j, k),
				nz(i, j, k),
				nx(i, j, k),
				alpha(i, j, k),
				m1,
				m2
			);
	}
	else
	{
		if (Q1 < 0.0)
		{
			vof1(i,j,k) = 
				calcV
				(
					nz(i, j, k),
					nx(i, j, k),
					ny(i, j, k),
					alpha(i,j,k),
					Q1,
					-Q1
				);
		}
		if (Q2 > 0.0) 
		{
			vof3(i,j,k) = 
				calcV
				(
					nz(i, j, k),
					nx(i, j, k),
					ny(i, j, k),
					alpha(i,j,k),
					1.0,
					Q2
				);
		}
		
		vof2(i,j,k) = 
			calcV
			(
				nz(i, j, k),
				nx(i, j, k),
				ny(i, j, k),
				alpha(i,j,k),
				m1,
				m2
			);
	}
}


double VOF_PLIC::calcV
(
    const double& m1,
    const double& m2,
    const double& m3,
    const double& alpha,
	double r0,
	double dr0
)
{
	//- Move origin to r0 along axis as given in Gueyffier (25)
	// if Q1 < 0: vof1 starts at r0=Q1 and goes to  0 --> dr0 = -Q1
	// if Q2 > 0: vof3 starts at r0= 1 and goes to Q2 --> dr0 = Q2
	// vof2 starts at r0=0 or r0=Q1 (if Q1 > 0) and ends at 1, 1-Q1 (if Q1 > 0), 1-Q1-abs(Q2) (if Q2 < 0)
	
	double al = alpha - m1*r0;

	
	//- Reflect parallelepiped such that m_i > 0
	
	al += max(0.0, -m1*dr0) + max(0.0, -m2) + max(0.0, -m3);
	
	
	//- Normalise plane equation: n1*x + n2*y + n3*z = al
	
	double tmp = fabs(m1)*dr0 + fabs(m2) + fabs(m3);
	if (tmp < 1e-10) return 0.0;	
	
	double n1 = fabs(m1)/(tmp + 1e-50);
	double n2 = fabs(m2)/(tmp + 1e-50);
	double n3 = fabs(m3)/(tmp + 1e-50);
	
	al = max(0.0, min(1.0, al/tmp));
	
	
	//- Limit alpha such that 0.0 < alpha < 0.5
	
	double al0 = min(al, 1.0 - al);
	
	
	//- Order b_i such that b1 < b2 < b3
	
	double b1 = min(n1*dr0, n2);
	double b3 = max(n1*dr0, n2);
	double b2 = n3;
	if (b2 < b1)
	{
		tmp = b1;
		b1 = b2;
		b2 = tmp;
	}
	else if (b2 > b3)
	{
		tmp = b3;
		b3 = b2;
		b2 = tmp;
	}
	double b12 = b1 + b2;
	double bm = min(b12, b3);
	
	
	//- Calculate new volume, assuming x0 = 0.0 and delta_x = 1.0 
	// Scardovelli p.233: alpha = al0, m1 = b1, m2 = m2, m3 = b3, m = bm, m12 = b12, m3 = b3, V = tmp
	
	double pr = max(6.0*b1*b2*b3, 1.0e-20);
	
	// if alpha < b1
	if (al0 < b1)
	{
		tmp = pow(al0, 3.0)/pr;
	}
	else if (al0 < b2)	// if b1 < alpha < b2
	{
		tmp = 0.5*al0*(al0 - b1)/(b2*b3 + 1e-20) + pow(b1, 3.0)/pr;
	}
	else if (al0 < bm)	// if b2 < alpha < bm
	{
		tmp = 
			(al0*al0*(3.0*b12 - al0) + b1*b1*(b1 - 3.0*al0) 
			+ b2*b2*(b2 - 3.0*al0))
			/pr;
	}
	else if (b12 < b3)	// if b12 < b3 --> bm = b12
	{
		tmp = (al0 - 0.5*bm)/(b3 + 1e-20);
	}
	else	// if b3 < b12 --> bm = b3
	{
		tmp = 
			(al0*al0*(3.0 - 2.0*al0) + b1*b1*(b1 - 3.0*al0) 
			+ b2*b2*(b2 - 3.0*al0) + b3*b3*(b3 - 3.0*al0))
			/pr;
	}
      
	double V = 0.0;  
	if (al <= 0.5)
	{
		V = tmp*dr0;
	}
	else
	{
		V = (1.0 - tmp)*dr0;
	}

	return V;
}


double VOF_PLIC::calcV2(lexer *p)
{
	double m1c1 = nx(i,j,k)*p->DXN[IP];
	double m2c2 = ny(i,j,k)*p->DYN[JP];
	double m3c3 = nz(i,j,k)*p->DZN[KP];
	
	double alpha_max = m1c1 + m2c2 + m3c3;
	
	if (fabs(alpha_max) > 1e-6)
	{
		double V = 
			1.0/(6.0*m1c1*m2c2*m3c3)*
			(
				pow(alpha(i,j,k),3.0) 
				- F3(alpha(i,j,k) - m1c1) 
				- F3(alpha(i,j,k) - m2c2) 
				- F3(alpha(i,j,k) - m3c3)
				+ F3(alpha(i,j,k) - alpha_max + m1c1) 
				+ F3(alpha(i,j,k) - alpha_max + m2c2) 
				+ F3(alpha(i,j,k) - alpha_max + m3c3)
			);

		return V;
	}
	else
	{
		return 0.0;
	}
}


double VOF_PLIC::F3(const double& x)
{	
	return (x <= 0.0) ? 0.0 : pow(x,3.0);
}
