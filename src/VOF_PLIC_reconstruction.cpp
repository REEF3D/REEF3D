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

void VOF_PLIC::reconstructPlane(fdm* a, lexer* p)
{
	//- Calculating interface normal vector n
	
	//calcNormalFO(a, p);
	//calcNormalWENO(a, p);
	//calcNormalPhi(a, p);
	calcNormalLS(a, p);


	//- Scale n_i according to cell size
	
	nx(i,j,k) *= p->DXN[IP];
	ny(i,j,k) *= p->DYN[JP];
	nz(i,j,k) *= p->DZN[KP];
	
	
    //- Ensure that n_i > 0
	
    int invx = -1;
    int invy = -1;
    int invz = -1;
    
    if (nx(i,j,k) < 0.0)
    {
        nx(i,j,k) *= -1.0;
        invx = 1;
    }
    if (ny(i,j,k) < 0.0)
    {
        ny(i,j,k) *= -1.0;
        invy = 1;
    }
    if (nz(i,j,k) < 0.0)
    {
        nz(i,j,k) *= -1.0;
        invz = 1;
    }


	//- Normalise plane

    double sum = nx(i,j,k) + ny(i,j,k) + nz(i,j,k);// + pow(p->DXN[IP]*p->DYN[JP]*p->DZN[KP],1.0/3.0)*1e-5;
	
	sum = sum > 1.0e-20 ? sum : 1.0e20;
	
    nx(i,j,k) /= sum;
    ny(i,j,k) /= sum;
    nz(i,j,k) /= sum;		


	//-  Calculating alpha from n and vof according to Scardovelli p.234
		
	alpha(i, j, k) = calcAlpha(a, nx(i, j, k), ny(i, j, k), nz(i, j, k));


	//- Return to original plane
	
	nx(i, j, k) *= -invx;
	ny(i, j, k) *= -invy;
	nz(i, j, k) *= -invz;
	
	alpha(i,j,k) += 
		min(0.0, nx(i, j, k)) + min(0.0, ny(i, j, k)) + min(0.0, nz(i, j, k));	
}


double VOF_PLIC::calcAlpha
(
	fdm* a,
	double& nx,
	double& ny,
	double& nz
)
{
	//- Assign n_i to m_i such that m1 < m2 < m3
	
	double tmp;
	double m1 = min(nx, ny);
	double m3 = max(nx, ny);
	double m2 = nz;
	
	if (m2 < m1)
	{
		tmp = m1;
		m1 = m2;
		m2 = tmp;
	}
	else if (m2 > m3) 
	{
		tmp = m3;
		m3 = m2;
		m2 = tmp;
	}
	

	//- Calculate ranges of functions V1, V2, V3
	
	double m12 = m1 + m2;
	double pr = max(6.0*m1*m2*m3, 1.0e-20);
	
	double V1  = pow(m1, 3.0)/pr;
	double V2  = V1 + 0.5*(m2 - m1)/(m3 + 1e-20);
	
	double V3, mm;
	if (m3 < m12)	// if m = m3
	{
		mm = m3;
		V3 = 
			(m3*m3*(3.0*m12 - m3) + m1*m1*(m1 - 3.0*m3) + m2*m2*(m2 - 3.0*m3))
			/pr;
	}
	else	// if m = m12
	{
		mm = m12;
		V3 = 0.5*mm/(m3 + 1e-20);
	}
    
	
	//- Limit vof such that 0 < vof < 0.5
	
	double vofLim = min(a->vof(i, j, k), 1.0 - a->vof(i, j, k));
     	 
	 
	//- Calculate alpha from vofLim and V_i

	double alpha;
	if (vofLim <= V1)	// 0 < V < V1
	{
		alpha = pow(pr*vofLim, 1.0/3.0);
	}
	else if (vofLim < V2)	// V1 < V < V2
	{
		alpha = 0.5*(m1 + sqrt(m1*m1 + 8.0*m2*m3*(vofLim-V1)));
	}
	else if (vofLim < V3)	// V2 < V < V3
	{
		double p = 2.0*m1*m2;
		double cs = cos(acos((1.5*m1*m2*(m12 - 2.0*m3*vofLim))/(p*sqrt(p)))/3.0);
		
		alpha = sqrt(p)*(sqrt(3.0*(1.0 - cs*cs)) - cs) + m12;
	}
	else if (m12 < m3)	// V32 < V < 0.5
	{
         alpha = m3*vofLim + 0.5*mm;
	}
	else	// V31 < V < 0.5
	{
         double p = m1*(m2 + m3) + m2*m3 - 0.25;
         double cs = cos(acos((1.5*m1*m2*m3*(0.5 - vofLim))/(p*sqrt(p)))/3.0);
		 
         alpha = sqrt(p)*(sqrt(3.0*(1.0 - cs*cs)) - cs) + 0.5;
	}


	//- Inverse result if 0.5 < vof < 1.0
	
	if (a->vof(i,j,k) > 0.5)  
	{
		alpha = 1.0 - alpha;
	}

	return alpha;
}
