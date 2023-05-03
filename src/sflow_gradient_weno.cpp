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

#include"sflow_gradient_weno.h"
#include"lexer.h"
#include"fdm2D.h"

sflow_gradient_weno::sflow_gradient_weno(lexer* p):tttw(13.0/12.0),fourth(1.0/4.0),third(1.0/3.0),
			sevsix(7.0/6.0),elvsix(11.0/6.0),sixth(1.0/6.0),fivsix(5.0/6.0),tenth(1.0/10.0),
			sixten(6.0/10.0),treten(3.0/10.0),epsilon(0.000001),smallnum(1.0e-20)
{
}

sflow_gradient_weno::~sflow_gradient_weno()
{
}

double sflow_gradient_weno::ddx(lexer *p,fdm2D *b, slice& f, int ipol, double advec)
{
    grad = 0.0;
	
	if(advec>0.0)
	{
	iqmin(p,b,f,ipol);
	is(f);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}
	

	if(advec<0.0)
	{
	iqmax(p,b,f,ipol);
	is(f);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

    //cout<<"GRAD: "<<grad<<"    "<<(f(i+1,j)-f(i,j))/p->DXM<<endl;
	return grad;
}

double sflow_gradient_weno::ddy(lexer *p,fdm2D *b, slice& f, int ipol, double advec)
{
    grad = 0.0;
	
	if(advec>0.0)
	{
	jqmin(p,b,f,ipol);
	is(f);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	
	}
	

	if(advec<0.0)
	{
	jqmax(p,b,f,ipol);
	is(f);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	
	}
	
	return grad;
}


void sflow_gradient_weno::iqmin(lexer *p,fdm2D *b, slice& f, int ipol)
{	
	q1 = (f(i-2,j) - f(i-3,j))/p->DXM;
	q2 = (f(i-1,j) - f(i-2,j))/p->DXM;
	q3 = (f(i,j)   - f(i-1,j))/p->DXM;
	q4 = (f(i+1,j) - f(i,j)  )/p->DXM;
	q5 = (f(i+2,j) - f(i+1,j))/p->DXM;
}

void sflow_gradient_weno::jqmin(lexer *p,fdm2D *b, slice& f, int ipol)
{
	q1 = (f(i,j-2) - f(i,j-3))/p->DXM;
	q2 = (f(i,j-1) - f(i,j-2))/p->DXM;
	q3 = (f(i,j)   - f(i,j-1))/p->DXM;
	q4 = (f(i,j+1) - f(i,j)  )/p->DXM;
	q5 = (f(i,j+2) - f(i,j+1))/p->DXM;
}

void sflow_gradient_weno::iqmax(lexer *p,fdm2D *b, slice& f, int ipol)
{
	q1 = (f(i+3,j) - f(i+2,j))/p->DXM;
	q2 = (f(i+2,j) - f(i+1,j))/p->DXM;
	q3 = (f(i+1,j) - f(i,j)  )/p->DXM;
	q4 = (f(i,j)   - f(i-1,j))/p->DXM;
	q5 = (f(i-1,j) - f(i-2,j))/p->DXM;
}

void sflow_gradient_weno::jqmax(lexer *p,fdm2D *b, slice& f, int ipol)
{
	q1 = (f(i,j+3) - f(i,j+2))/p->DXM;
	q2 = (f(i,j+2) - f(i,j+1))/p->DXM;
	q3 = (f(i,j+1) - f(i,j)  )/p->DXM;
	q4 = (f(i,j)   - f(i,j-1))/p->DXM;
	q5 = (f(i,j-1) - f(i,j-2))/p->DXM;
}

void sflow_gradient_weno::is(slice& f)
{
	is1 = tttw*pow(q1 - 2.0*q2 + q3, 2.0) + fourth*pow(q1 - 4.0*q2 + 3.0*q3, 2.0);
	is2 = tttw*pow(q2 - 2.0*q3 + q4, 2.0) + fourth*pow(q2 - q4, 2.0);
	is3 = tttw*pow(q3 - 2.0*q4 + q5, 2.0) + fourth*pow(3.0*q3 - 4.0*q4 + q5, 2.0);
}

void sflow_gradient_weno::alpha()
{
	alpha1=tenth/pow(epsilon+is1,2.0);
	alpha2=sixten/pow(epsilon+is2,2.0);
	alpha3=treten/pow(epsilon+is3,2.0);
}

void sflow_gradient_weno::weight()
{
	w1=alpha1/(alpha1+alpha2+alpha3);
	w2=alpha2/(alpha1+alpha2+alpha3);
	w3=alpha3/(alpha1+alpha2+alpha3);
}
