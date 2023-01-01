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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"limo3.h"
#include"lexer.h"
#include"fdm.h"

limo3::limo3 (lexer *p) : delta(p->DXM), radius (0.1), eps(1.0e-9*p->DXM)
{
}

limo3::~limo3()
{

}

double limo3::iphi(field& b,int n1, int n2, int q1, int q2)
{	
	d1 = b(i+n1,j,k)-b(i+n2,j,k);
	d2 = b(i+q1,j,k)-b(i+q2,j,k);
	
    r=d1/(fabs(d2)>1.0e-10?d2:1.0e20);
	
	eta = (d1*d1 + d2*d2)/pow(radius*delta,2.0);
	
	phihat = max(0.0, min((2.0+r)/3.0, max(-0.5*r, min(2.0*r, (2.0+r)/3.0, 1.6))));
	
	
	if(eta<=1.0-eps)
	phi = (2.0+r)/3.0;
	
	if(eta>=1.0+eps)
	phi = phihat;
	
	if(eta>1.0-eps && eta<1.0+eps)
	phi = 0.5 * ((1.0 - (eta-1.0)/eps)*(2.0+r)/3.0  + (1.0 + (eta-1.0)/eps)*phihat);

    return phi;
}

double limo3::jphi(field& b,int n1, int n2, int q1, int q2)
{
	d1=b(i,j+n1,k)-b(i,j+n2,k);
    d2=b(i,j+q1,k)-b(i,j+q2,k);
    
	r=d1/(fabs(d2)>1.0e-10?d2:1.0e20);
	
	eta = (d1*d1 + d2*d2)/pow(radius*delta,2.0);
	
	phihat = max(0.0, min((2.0+r)/3.0, max(-0.5*r, min(2.0*r, (2.0+r)/3.0, 1.6))));
	
	
	if(eta<=1.0-eps)
	phi = (2.0+r)/3.0;
	
	if(eta>=1.0+eps)
	phi = phihat;
	
	if(eta>1.0-eps && eta<1.0+eps)
	phi = 0.5 * ((1.0 - (eta-1.0)/eps)*(2.0+r)/3.0  + (1.0 + (eta-1.0)/eps)*phihat);

    return phi;
}

double limo3::kphi(field& b,int n1, int n2, int q1, int q2)
{		
	d1=b(i,j,k+n1)-b(i,j,k+n2);
    d2=b(i,j,k+q1)-b(i,j,k+q2);
	
    r=d1/(fabs(d2)>1.0e-10?d2:1.0e20);
	
	eta = (d1*d1 + d2*d2)/pow(radius*delta,2.0);
	
	phihat = max(0.0, min((2.0+r)/3.0, max(-0.5*r, min(2.0*r, (2.0+r)/3.0, 1.6))));
	
	
	if(eta<=1.0-eps)
	phi = (2.0+r)/3.0;
	
	if(eta>=1.0+eps)
	phi = phihat;
	
	if(eta>1.0-eps && eta<1.0+eps)
	phi = 0.5 * ((1.0 - (eta-1.0)/eps)*(2.0+r)/3.0  + (1.0 + (eta-1.0)/eps)*phihat);

    return phi;
}


double limo3::min(double val1,double val2,double val3)
{
	double mini;

	mini=val1;

	if(mini>val2)
	mini=val2;

	if(mini>val3)
	mini=val3;

	if(mini<0.0)
	mini=0.0;

	return mini;
}

double limo3::min(double val1,double val2)
{
	double mini;

	mini=val1;

	if(mini>val2)
	mini=val2;

	if(mini<0.0)
	mini=0.0;

	return mini;
}

double limo3::max(double val1,double val2,double val3)
{
	double maxi;

	maxi=val1;

	if(maxi<val2)
	maxi=val2;

	if(maxi<val3)
	maxi=val3;

	if(maxi<0.0)
	maxi=0.0;

	return maxi;
}

double limo3::max(double val1,double val2)
{
	double maxi;

	maxi=val1;

	if(maxi<val2)
	maxi=val2;

	if(maxi<0.0)
	maxi=0.0;

	return maxi;
}

