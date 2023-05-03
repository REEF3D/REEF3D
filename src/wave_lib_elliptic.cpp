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

#include"wave_lib_elliptic.h"
#include"lexer.h"
#include"ghostcell.h"

wave_lib_elliptic::wave_lib_elliptic() : epsi(1.0e-19)
{ 
}

wave_lib_elliptic::~wave_lib_elliptic()
{
}


void wave_lib_elliptic::elliptic(lexer *p, double u, double &sn, double &cn, double &dn)
{ 
	double res = 1.0e-10;
	double sinu,cosu,r; 
	int maxiter=18;
	
	if(modulus<epsi)
	{
	sn = sin(u);
	cn = cos(u);
	dn = 1.0;
	}
	
	else if(1.0-modulus<epsi)
	{
	sn = tanh(u);
	cn = 1.0/cosh(u);
	dn = cn;
	}
	
	else
	{
	double mu[20], nu[20], c[20], d[20];
	
	mu[0] = 1.0;
	nu[0] = sqrt(1.0 - modulus);
		
	int	n=0;
		// first reecursion
		while(fabs((mu[n] - nu[n])/(mu[n] + nu[n] + 1.0e-20)) > res)  // change
		{
			mu[n+1] = 0.5*(mu[n] + nu[n]);
			nu[n+1] = sqrt(mu[n]*nu[n]);
		
		++n;
		
		if(n>=maxiter)
		break;
		}	
		
	sinu = sin(u*mu[n]);
	cosu = cos(u*mu[n]);
	
		// second recursion
		if(fabs(sinu) < fabs(cosu))
		{ 
		c[n] = mu[n] * (sinu/cosu);
		d[n] = 1.0;
		
			while(n>0)
			{
			--n;
			
			c[n] = d[n+1] * c[n+1];
			r = (c[n+1]*c[n+1])/mu[n+1];
			
			d[n] = (r + nu[n])/(r + mu[n]); 
			}
			
		dn = sqrt(1.0 - modulus)/d[n];
		cn = dn*SIGN(cosu)/sqrt(1.0 + c[n]*c[n]);
		sn = cn*(c[n]/sqrt(1.0-modulus));
		}
		
		else
		{
		c[n] = mu[n] * (cosu/sinu);
		d[n] = 1.0;
		
			while(n>0)
			{
			--n;
			
			c[n] = d[n+1]*c[n+1];
			r = (c[n+1]*c[n+1])/mu[n+1];
			
			d[n] = (r + nu[n])/(r + mu[n]);				
			}
			
		dn = d[n];
		sn = SIGN(sinu)/sqrt(1.0 + c[n]*c[n]);
		cn = c[n]*sn;		
		}
		
	}
    
}

double wave_lib_elliptic::K_elliptic_1(double m)
{
	// Abramowitz, 17.3.33
	
	double a[3],b[3];
	double m1 = 1.0 - m;
	double K;
	
	
	a[0] = 1.3862944;
	a[1] = 0.1119723;
	a[2] = 0.0725296;

	b[0] = 0.5;
	b[1] = 0.1212478;
	b[2] = 0.0288729;
	
	K = a[0] + a[1]*m1 + a[2]*m1*m1 
	  +(b[0] + b[1]*m1 + b[2]*m1*m1)*log(1.0/m1); 
	  
	return K; 	
}

double wave_lib_elliptic::E_elliptic_1(double m)
{
	// Abramowitz, 17.3.35
	
	double a[3],b[3];
	double m1 = 1.0 - m;
	double E;	
	
	a[0] = 1.0;
	a[1] = 0.4630151;
	a[2] = 0.1077812;
	
	b[0] = 0.0;
	b[1] = 0.2452727;
	b[2] = 0.0412496;
	
	E = a[0] + a[1]*m1 + a[2]*m1*m1
	  +(b[0] + b[1]*m1 + b[2]*m1*m1)*log(1.0/m1); 
 	
	return E;
}

double wave_lib_elliptic::K_elliptic_5(double m)
{
	// Abramowitz, 17.3.33
	
	double a[3],b[3];
	double m1 = 1.0 - m;
	double K;
	
	
	a[0] = 1.3862944;
	a[1] = 0.1119723;
	a[2] = 0.0725296;

	b[0] = 0.5;
	b[1] = 0.1212478;
	b[2] = 0.0288729;
	
	K = a[0] + a[1]*m1 + a[2]*m1*m1 
	  +(b[0] + b[1]*m1 + b[2]*m1*m1)*log(1.0/m1); 
	  
	return K; 	
}

double wave_lib_elliptic::E_elliptic_5(double m)
{
	// Abramowitz, 17.3.35
	
	double a[3],b[3];
	double m1 = 1.0 - m;
	double E;	
	
	a[0] = 1.0;
	a[1] = 0.4630151;
	a[2] = 0.1077812;
	
	b[0] = 0.0;
	b[1] = 0.2452727;
	b[2] = 0.0412496;
	
	E = a[0] + a[1]*m1 + a[2]*m1*m1
	  +(b[0] + b[1]*m1 + b[2]*m1*m1)*log(1.0/m1); 
 	
	return E;
}

double wave_lib_elliptic::K_elliptic(double m)
{
	// Abramowitz, 17.3.34
	
	double a[5],b[5];
	double m1 = 1.0 - m;
	double K;
	
	
	a[0] = 1.38629;
	a[1] = 0.09666;
	a[2] = 0.03590;
	a[3] = 0.03742;
	a[4] = 0.01451;
	
	b[0] = 0.5;
	b[1] = 0.12498;
	b[2] = 0.06880;
	b[3] = 0.03328;
	b[4] = 0.00441;
	
	K = a[0] + a[1]*m1 + a[2]*m1*m1 + a[3]*m1*m1*m1 + a[4]*m1*m1*m1*m1
	  +(b[0] + b[1]*m1 + b[2]*m1*m1 + b[3]*m1*m1*m1 + b[4]*m1*m1*m1*m1)*log(1.0/m1); 
	  
	return K; 	
}

double wave_lib_elliptic::E_elliptic(double m)
{
	// Abramowitz, 17.3.36
	
	double a[5],b[5];
	double m1 = 1.0 - m;
	double E;	
	
	a[0] = 1.0;
	a[1] = 0.44325;
	a[2] = 0.06260;
	a[3] = 0.04757;
	a[4] = 0.01736;
	
	b[0] = 0.0;
	b[1] = 0.24998;
	b[2] = 0.09200;
	b[3] = 0.04069;
	b[4] = 0.00526;
	
	E = a[0] + a[1]*m1 + a[2]*m1*m1 + a[3]*m1*m1*m1 + a[4]*m1*m1*m1*m1
	  +(b[0] + b[1]*m1 + b[2]*m1*m1 + b[3]*m1*m1*m1 + b[4]*m1*m1*m1*m1)*log(1.0/m1); 
 	
	return E;
}
