/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"air_entrainment_ibc.h"
#include"lexer.h"
#include"fdm.h"
#include"turbulence.h"

air_entrainment_ibc::air_entrainment_ibc(lexer *p)
{

}

air_entrainment_ibc::~air_entrainment_ibc()
{
}

void air_entrainment_ibc::air_entrainment_ibc_start(lexer* p,fdm* a,ghostcell *pgc,turbulence *pturb,field& conc)
{
	double epsi=3.1*p->dx;
	double dirac,area;
	double dstx,dsty,dstz,dnorm;
	double x_norm,y_norm,z_norm;
	double Pd,Pt,Lt;
	double Ca=0.5;
	double vel;
	int count;
	
	count=0;
	LOOP
	{	
		if(a->phi(i,j,k)<epsi && a->phi(i,j,k)>-0.6*p->dx)
		{
		dirac = (0.5/epsi)*(1.0 + cos((PI*a->phi(i,j,k))/epsi));
		
		
		dstx = (a->phi(i+1,j,k)-a->phi(i-1,j,k))/(2.0*p->dx);
		dsty = (a->phi(i,j+1,k)-a->phi(i,j-1,k))/(2.0*p->dx);
		dstz = (a->phi(i,j,k+1)-a->phi(i,j,k-1))/(2.0*p->dx);
		
		dnorm=sqrt(dstx*dstx + dsty*dsty + dstz*dstz);
		x_norm = dstx/(dnorm>1.0e-20?dnorm:1.0e20);
		y_norm = dsty/(dnorm>1.0e-20?dnorm:1.0e20);
		z_norm = dstz/(dnorm>1.0e-20?dnorm:1.0e20);
		
		vel = sqrt(a->u(i,j,k)*a->u(i,j,k) + a->v(i,j,k)*a->v(i,j,k) + a->w(i,j,k)*a->w(i,j,k));
		
		area =  pow(p->dx,3.0) * dirac *dnorm;
		
		Lt = sqrt(1.5)*pow(fabs(pturb->kinval(i,j,k)),0.5)/((pturb->epsval(i,j,k))>(1.0e-20)?(pturb->epsval(i,j,k)):(1.0e20));
		
		Pd = z_norm*Lt;
		Pt = pturb->kinval(i,j,k);
		
			if(Pt>Pd && vel>0.8)
			{ 
				a->rhsvec.V[count] += dirac*Ca*area*sqrt(2.0*(Pt-Pd)) * (1.0/pow(p->dx,2.0)) ;
			}
		
		
		}
		
	++count;	
	}
}

