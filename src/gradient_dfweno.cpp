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

#include"gradient.h"
#include"fdm.h"
#include<math.h>

double gradient::dfwenox(fdm* a, field& b, double uw)
{				
		i-=1;
		fu1 = fx(a,b,uw);
		i+=1;
		
		fu2 = fx(a,b,uw);
		
		ddx = (fu2-fu1)/dx;
		
		return ddx;
}

double gradient::dfwenoy(fdm* a, field& b, double uw)
{				
		j-=1;
		fv1 = fy(a,b,uw);
		j+=1;
		
		fv2 = fy(a,b,uw);
		
		ddy = (fv2-fv1)/dx;
		
		return ddy;
}

double gradient::dfwenoz(fdm* a, field& b, double uw)
{				
		k-=1;
		fw1 = fz(a,b,uw);
		k+=1;
		
		fw2 = fz(a,b,uw);
		
		ddz = (fw2-fw1)/dx;
		
		return ddz;
}
		
double gradient::fx(fdm* a, field& b, double uw)
{
    grad = 0.0;

	if(uw>0.0)
	{
	iqfmin(b);
	isf(b);
	alphaf();
	weightf();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	if(uw<0.0)
	{
	iqfmax(b);
	isf(b);
	alphaf();
	weightf();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	return grad;
}

double gradient::fy(fdm* a, field& b, double uw)
{
    grad = 0.0;

	if(uw>0.0)
	{
	jqfmin(b);
	isf(b);
	alphaf();
	weightf();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	if(uw<0.0)
	{
	jqfmax(b);
	isf(b);
	alphaf();
	weightf();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}
	
	return grad;
}

double gradient::fz(fdm* a, field& b, double uw)
{
    grad = 0.0;

	if(uw>0.0)
	{
	kqfmin(b);
	isf(b);
	alphaf();
	weightf();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	if(uw<0.0)
	{
	kqfmax(b);
	isf(b);
	alphaf();
	weightf();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	return grad;
}

void gradient::iqfmin(field& f)
{	
	q1 = f(i-2,j,k);
	q2 = f(i-1,j,k);
	q3 = f(i,j,k);
	q4 = f(i+1,j,k);
	q5 = f(i+2,j,k);
}

void gradient::jqfmin(field& f)
{
	q1 = f(i,j-2,k);
	q2 = f(i,j-1,k);
	q3 = f(i,j,k);
	q4 = f(i,j+1,k);
	q5 = f(i,j+2,k);
}

void gradient::kqfmin(field& f)
{
	q1 = f(i,j,k-2);
	q2 = f(i,j,k-1);
	q3 = f(i,j,k);
	q4 = f(i,j,k+1);
	q5 = f(i,j,k+2);
}

void gradient::iqfmax(field& f)
{
	q1 = f(i+3,j,k);
	q2 = f(i+2,j,k);
	q3 = f(i+1,j,k);
	q4 = f(i,j,k);
	q5 = f(i-1,j,k);
}

void gradient::jqfmax(field& f)
{
	q1 = f(i,j+3,k);
	q2 = f(i,j+2,k);
	q3 = f(i,j+1,k);
	q4 = f(i,j,k);
	q5 = f(i,j-1,k);
}

void gradient::kqfmax(field& f)
{
	q1 = f(i,j,k+3);
	q2 = f(i,j,k+2);
	q3 = f(i,j,k+1);
	q4 = f(i,j,k);
	q5 = f(i,j,k-1);
}

void gradient::isf(field& f)
{
	is1 = tttw*pow(q1-2.0*q2+q3, 2.0) + fourth*pow(q1-4.0*q2+3.0*q3, 2.0);
	is2 = tttw*pow(q2-2.0*q3+q4, 2.0) + fourth*pow(q2-q4, 2.0);
	is3 = tttw*pow(q3-2.0*q4+q5, 2.0) + fourth*pow(3.0*q3-4.0*q4+q5, 2.0);
}

void gradient::alphaf()
{
	alpha1=tenth*pow(epsilon+is1,-2.0);
	alpha2=sixten*pow(epsilon+is2,-2.0);
	alpha3=treten*pow(epsilon+is3,-2.0);
}

void gradient::weightf()
{
	w1=alpha1/(alpha1+alpha2+alpha3);
	w2=alpha2/(alpha1+alpha2+alpha3);
	w3=alpha3/(alpha1+alpha2+alpha3);
}
