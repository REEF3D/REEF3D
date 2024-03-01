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

double nhflow_gradient::dwenox(double *F, double uw)
{
	grad=0.0;

	if(uw>=0.0)
	{
	iqmin(F,p->DXD);
	is();
	alpha();
	weight();
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	if(uw<0.0)
	{
	iqmax(F,p->DXD);
	is();
	alpha();
	weight();
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	return grad;
}

double nhflow_gradient::dwenoy(double *F, double uw)
{
	grad=0.0;

	if(uw>=0.0)
	{
	jqmin(F,p->DYD);
	is();
	alpha();
	weight();
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	if(uw<0.0)
	{
	jqmax(F,p->DYD);
	is();
	alpha();
	weight();
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	return grad;
}

void nhflow_gradient::iqmin(double *F, double delta)
{
    q1 = (F[Im2JK] - F[Im3JK])/delta;
	q2 = (F[Im1JK] - F[Im2JK])/delta;
	q3 = (F[IJK]  -  F[Im1JK])/delta;
	q4 = (F[Ip1JK] - F[IJK]  )/delta;
	q5 = (F[Ip2JK] - F[Ip1JK])/delta;
}

void nhflow_gradient::jqmin(double *F, double delta)
{
	q1 = (F[IJm2K] - F[IJm3K])/delta;
	q2 = (F[IJm1K] - F[IJm2K])/delta;
	q3 = (F[IJK]  -  F[IJm1K])/delta;
	q4 = (F[IJp1K] - F[IJK]  )/delta;
	q5 = (F[IJp2K] - F[IJp1K])/delta;
}

void nhflow_gradient::iqmax(double *F, double delta)
{
    q1 = (F[Ip3JK]  - F[Ip2JK])/delta;
	q2 = (F[Ip2JK]  - F[Ip1JK])/delta;
	q3 = (F[Ip1JK]  - F[IJK])/delta;
	q4 = (F[IJK]    - F[Im1JK])/delta;
	q5 = (F[Im1JK]  - F[Im2JK])/delta;
}

void nhflow_gradient::jqmax(double *F, double delta)
{
	q1 = (F[IJp3K]  - F[IJp2K])/delta;
	q2 = (F[IJp2K]  - F[IJp2K])/delta;
	q3 = (F[IJp1K]  - F[IJK])/delta;
	q4 = (F[IJK]    - F[IJm1K])/delta;
	q5 = (F[IJm1K]  - F[IJm2K])/delta;
}


