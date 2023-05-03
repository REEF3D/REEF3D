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
#include"slice.h"

double gradient::dslwenox(fdm* a, slice& b, double uw)
{
	grad=0.0;

	if(uw>0.0)
	{
	iqminsl(b,dx);
	is();
	alpha();
	weight();
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	if(uw<0.0)
	{
	iqmaxsl(b,dx);
	is();
	alpha();
	weight();
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	return grad;
}

double gradient::dslwenoy(fdm* a, slice& b, double uw)
{
	grad=0.0;

	if(uw>0.0)
	{
	jqminsl(b,dx);
	is();
	alpha();
	weight();
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	if(uw<0.0)
	{
	jqmaxsl(b,dx);
	is();
	alpha();
	weight();
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	return grad;
}

void gradient::iqminsl(slice& f, double delta)
{
	q1 = (f(i-2,j) - f(i-3,j))/delta;
	q2 = (f(i-1,j) - f(i-2,j))/delta;
	q3 = (f(i,j)   - f(i-1,j))/delta;
	q4 = (f(i+1,j) - f(i,j)  )/delta;
	q5 = (f(i+2,j) - f(i+1,j))/delta;
}

void gradient::jqminsl(slice& f, double delta)
{
	q1 = (f(i,j-2) - f(i,j-3))/delta;
	q2 = (f(i,j-1) - f(i,j-2))/delta;
	q3 = (f(i,j)   - f(i,j-1))/delta;
	q4 = (f(i,j+1) - f(i,j)  )/delta;
	q5 = (f(i,j+2) - f(i,j+1))/delta;
}

void gradient::iqmaxsl(slice& f, double delta)
{
	q1 = (f(i+3,j) - f(i+2,j))/delta;
	q2 = (f(i+2,j) - f(i+1,j))/delta;
	q3 = (f(i+1,j) - f(i,j)  )/delta;
	q4 = (f(i,j)   - f(i-1,j))/delta;
	q5 = (f(i-1,j) - f(i-2,j))/delta;
}

void gradient::jqmaxsl(slice& f, double delta)
{
	q1 = (f(i,j+3) - f(i,j+2))/delta;
	q2 = (f(i,j+2) - f(i,j+1))/delta;
	q3 = (f(i,j+1) - f(i,j)  )/delta;
	q4 = (f(i,j)   - f(i,j-1))/delta;
	q5 = (f(i,j-1) - f(i,j-2))/delta;
}
