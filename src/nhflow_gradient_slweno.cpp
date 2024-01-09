/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"nhflow_gradient.h"
#include"lexer.h"
#include"slice.h"
#include<math.h>

double nhflow_gradient::dslwenox(slice& b, double uw)
{
	grad=0.0;

	if(uw>=0.0)
	{
	iqminsl(b,p->DXD);
	is();
	alpha();
	weight();
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	if(uw<0.0)
	{
	iqmaxsl(b,p->DXD);
	is();
	alpha();
	weight();
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	return grad;
}

double nhflow_gradient::dslwenoy(slice& b, double uw)
{
	grad=0.0;

	if(uw>=0.0)
	{
	jqminsl(b,p->DYD);
	is();
	alpha();
	weight();
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	if(uw<0.0)
	{
	jqmaxsl(b,p->DYD);
	is();
	alpha();
	weight();
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	return grad;
}

void nhflow_gradient::iqminsl(slice& f, double delta)
{
	q1 = (f(i-2,j) - f(i-3,j))/delta;
	q2 = (f(i-1,j) - f(i-2,j))/delta;
	q3 = (f(i,j)   - f(i-1,j))/delta;
	q4 = (f(i+1,j) - f(i,j)  )/delta;
	q5 = (f(i+2,j) - f(i+1,j))/delta;
}

void nhflow_gradient::jqminsl(slice& f, double delta)
{
	q1 = (f(i,j-2) - f(i,j-3))/delta;
	q2 = (f(i,j-1) - f(i,j-2))/delta;
	q3 = (f(i,j)   - f(i,j-1))/delta;
	q4 = (f(i,j+1) - f(i,j)  )/delta;
	q5 = (f(i,j+2) - f(i,j+1))/delta;
}

void nhflow_gradient::iqmaxsl(slice& f, double delta)
{
	q1 = (f(i+3,j) - f(i+2,j))/delta;
	q2 = (f(i+2,j) - f(i+1,j))/delta;
	q3 = (f(i+1,j) - f(i,j)  )/delta;
	q4 = (f(i,j)   - f(i-1,j))/delta;
	q5 = (f(i-1,j) - f(i-2,j))/delta;
}

void nhflow_gradient::jqmaxsl(slice& f, double delta)
{
	q1 = (f(i,j+3) - f(i,j+2))/delta;
	q2 = (f(i,j+2) - f(i,j+1))/delta;
	q3 = (f(i,j+1) - f(i,j)  )/delta;
	q4 = (f(i,j)   - f(i,j-1))/delta;
	q5 = (f(i,j-1) - f(i,j-2))/delta;
}

void nhflow_gradient::is()
{
	is1 = tttw*pow(q1-2.0*q2+q3,2.0) + fourth*pow((q1-4.0*q2+3.0*q3), 2.0);
	is2 = tttw*pow(q2-2.0*q3+q4, 2.0) + fourth*pow((q2-q4), 2.0);
	is3 = tttw*pow(q3-2.0*q4+q5, 2.0) + fourth*pow((3.0*q3-4.0*q4+q5), 2.0);
}

void nhflow_gradient::alpha()
{
	alpha1=tenth*pow(epsilon+is1,-2.0);
	alpha2=sixten*pow(epsilon+is2,-2.0);
	alpha3=treten*pow(epsilon+is3,-2.0);
}

void nhflow_gradient::weight()
{
	w1=alpha1/(alpha1+alpha2+alpha3);
	w2=alpha2/(alpha1+alpha2+alpha3);
	w3=alpha3/(alpha1+alpha2+alpha3);
}
