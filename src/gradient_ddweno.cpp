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
#include"lexer.h"
#include"fdm.h"
#include<math.h>

double gradient::ddwenox(fdm* a, field& b, double uw)
{
	grad=0.0;

	if(uw>0.0)
	{
	iqmin(b,dx);
	is();
	alpha();
	weight();
 //   if(p->mpirank==0)
//cout<<i<<" "<<j<<" "<<k<<" weigths: "<<w3<<" "<<w2<<" "<<w1<<endl;


	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}


	if(uw<0.0)
	{
	iqmax(b,dx);
	is();
	alpha();
	weight();
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	return grad;
}

double gradient::ddwenoy(fdm* a, field& b, double uw)
{
	grad=0.0;

	if(uw>0.0)
	{
	jqmin(b,dx);
	is();
	alpha();
	weight();
    
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	if(uw<0.0)
	{
	jqmax(b,dx);
	is();
	alpha();
	weight();
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	return grad;

}

double gradient::ddwenoz(fdm* a, field& b, double uw)
{
	grad=0.0;

	if(uw>0.0)
	{
	kqmin(b,dx);
	is();
	alpha();
	weight();
    
//        if(p->mpirank==0)
//cout<<i<<" "<<j<<" "<<k<<" weigths: "<<is3<<" "<<is2<<" "<<is1<<" weigths: "<<w3<<" "<<w2<<" "<<w1<<endl;

	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}


	if(uw<0.0)
	{
	kqmax(b,dx);
	is();
	alpha();
	weight();
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}
	return grad;
}

void gradient::iqmin(field& f, double delta)
{
	q1 = (f(i-2,j,k) - f(i-3,j,k))/delta;
	q2 = (f(i-1,j,k) - f(i-2,j,k))/delta;
	q3 = (f(i,j,k)   - f(i-1,j,k))/delta;
	q4 = (f(i+1,j,k) - f(i,j,k)  )/delta;
	q5 = (f(i+2,j,k) - f(i+1,j,k))/delta;
}

void gradient::jqmin(field& f, double delta)
{
	q1 = (f(i,j-2,k) - f(i,j-3,k))/delta;
	q2 = (f(i,j-1,k) - f(i,j-2,k))/delta;
	q3 = (f(i,j,k)   - f(i,j-1,k))/delta;
	q4 = (f(i,j+1,k) - f(i,j,k)  )/delta;
	q5 = (f(i,j+2,k) - f(i,j+1,k))/delta;
}

void gradient::kqmin(field& f, double delta)
{

	q1 = (f(i,j,k-2) - f(i,j,k-3))/delta;
	q2 = (f(i,j,k-1) - f(i,j,k-2))/delta;
	q3 = (f(i,j,k)   - f(i,j,k-1))/delta;
	q4 = (f(i,j,k+1) - f(i,j,k)  )/delta;
	q5 = (f(i,j,k+2) - f(i,j,k+1))/delta;

}

void gradient::iqmax(field& f, double delta)
{
	q1 = (f(i+3,j,k) - f(i+2,j,k))/delta;
	q2 = (f(i+2,j,k) - f(i+1,j,k))/delta;
	q3 = (f(i+1,j,k) - f(i,j,k)  )/delta;
	q4 = (f(i,j,k)   - f(i-1,j,k))/delta;
	q5 = (f(i-1,j,k) - f(i-2,j,k))/delta;
}

void gradient::jqmax(field& f, double delta)
{
	q1 = (f(i,j+3,k) - f(i,j+2,k))/delta;
	q2 = (f(i,j+2,k) - f(i,j+1,k))/delta;
	q3 = (f(i,j+1,k) - f(i,j,k)  )/delta;
	q4 = (f(i,j,k)   - f(i,j-1,k))/delta;
	q5 = (f(i,j-1,k) - f(i,j-2,k))/delta;
}

void gradient::kqmax(field& f, double delta)
{
	q1 = (f(i,j,k+3) - f(i,j,k+2))/delta;
	q2 = (f(i,j,k+2) - f(i,j,k+1))/delta;
	q3 = (f(i,j,k+1) - f(i,j,k)  )/delta;
	q4 = (f(i,j,k)   - f(i,j,k-1))/delta;
	q5 = (f(i,j,k-1) - f(i,j,k-2))/delta;
}

void gradient::is()
{
	is1 = tttw*pow(q1-2.0*q2+q3,2.0) + fourth*pow((q1-4.0*q2+3.0*q3), 2.0);
	is2 = tttw*pow(q2-2.0*q3+q4, 2.0) + fourth*pow((q2-q4), 2.0);
	is3 = tttw*pow(q3-2.0*q4+q5, 2.0) + fourth*pow((3.0*q3-4.0*q4+q5), 2.0);
}

void gradient::alpha()
{
	alpha1=tenth*pow(epsilon+is1,-2.0);
	alpha2=sixten*pow(epsilon+is2,-2.0);
	alpha3=treten*pow(epsilon+is3,-2.0);
}

void gradient::weight()
{
	w1=alpha1/(alpha1+alpha2+alpha3);
	w2=alpha2/(alpha1+alpha2+alpha3);
	w3=alpha3/(alpha1+alpha2+alpha3);
}
