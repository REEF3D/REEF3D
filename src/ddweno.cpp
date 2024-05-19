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

#include"ddweno.h"
#include"fdm.h"
#include"lexer.h"
#include"ghostcell.h"
#include"vec.h"
#include"cpt.h"

ddweno::ddweno(lexer* p):tttw(13.0/12.0),fourth(1.0/4.0),third(1.0/3.0),
			sevsix(7.0/6.0),elvsix(11.0/6.0),sixth(1.0/6.0),fivsix(5.0/6.0),tenth(1.0/10.0),
			sixten(6.0/10.0),treten(3.0/10.0),epsilon(0.000001),smallnum(1.0e-20),dx(p->DXM)
			
{
}

ddweno::~ddweno()
{

}

double ddweno::ddwenox(fdm* a, vec& b, double uw, cpt &C)
{

	grad=0.0;

	if(uw>0.0)
	{
	iqmin(b,dx,C);
	is();
	alpha();
	weight();
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	if(uw<0.0)
	{
	iqmax(b,dx,C);
	is();
	alpha();
	weight();
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	return grad;
}

double ddweno::ddwenoy(fdm* a, vec& b, double uw, cpt &C)
{
	grad=0.0;

	if(uw>0.0)
	{
	jqmin(b,dx,C);
	is();
	alpha();
	weight();
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	if(uw<0.0)
	{
	jqmax(b,dx,C);
	is();
	alpha();
	weight();
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	return grad;
}

double ddweno::ddwenoz(fdm* a, vec& b, double uw, cpt &C)
{
	grad=0.0;

	if(uw>0.0)
	{
	kqmin(b,dx,C);
	is();
	alpha();
	weight();
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}


	if(uw<0.0)
	{
	kqmax(b,dx,C);
	is();
	alpha();
	weight();
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}
	return grad;
}

void ddweno::iqmin(vec& f, double delta, cpt &C)
{
	q1 = (f.V[Im2_J_K] - f.V[Im3_J_K])/delta;
	q2 = (f.V[Im1_J_K] - f.V[Im2_J_K])/delta;
	q3 = (f.V[I_J_K]   - f.V[Im1_J_K])/delta;
	q4 = (f.V[Ip1_J_K] - f.V[I_J_K]  )/delta;
	q5 = (f.V[Ip2_J_K] - f.V[Ip1_J_K])/delta;
}

void ddweno::jqmin(vec& f, double delta, cpt &C)
{
	q1 = (f.V[I_Jm2_K] - f.V[I_Jm3_K])/delta;
	q2 = (f.V[I_Jm1_K] - f.V[I_Jm2_K])/delta;
	q3 = (f.V[I_J_K]   - f.V[I_Jm1_K])/delta;
	q4 = (f.V[I_Jp1_K] - f.V[I_J_K]  )/delta;
	q5 = (f.V[I_Jp2_K] - f.V[I_Jp1_K])/delta;
}

void ddweno::kqmin(vec& f, double delta, cpt &C)
{
	q1 = (f.V[I_J_Km2] - f.V[I_J_Km3])/delta;
	q2 = (f.V[I_J_Km1] - f.V[I_J_Km2])/delta;
	q3 = (f.V[I_J_K]   - f.V[I_J_Km1])/delta;
	q4 = (f.V[I_J_Kp1] - f.V[I_J_K]  )/delta;
	q5 = (f.V[I_J_Kp2] - f.V[I_J_Kp1])/delta;
}

void ddweno::iqmax(vec& f, double delta, cpt &C)
{
	q1 = (f.V[Ip3_J_K] - f.V[Ip2_J_K])/delta;
	q2 = (f.V[Ip2_J_K] - f.V[Ip1_J_K])/delta;
	q3 = (f.V[Ip1_J_K] - f.V[I_J_K]  )/delta;
	q4 = (f.V[I_J_K]   - f.V[Im1_J_K])/delta;
	q5 = (f.V[Im1_J_K] - f.V[Im2_J_K])/delta;
}

void ddweno::jqmax(vec& f, double delta, cpt &C)
{
	q1 = (f.V[I_Jp3_K] - f.V[I_Jp2_K])/delta;
	q2 = (f.V[I_Jp2_K] - f.V[I_Jp1_K])/delta;
	q3 = (f.V[I_Jp1_K] - f.V[I_J_K]  )/delta;
	q4 = (f.V[I_J_K]   - f.V[I_Jm1_K])/delta;
	q5 = (f.V[I_Jm1_K] - f.V[I_Jm2_K])/delta;
}

void ddweno::kqmax(vec& f, double delta, cpt &C)
{
	q1 = (f.V[I_J_Kp3] - f.V[I_J_Kp2])/delta;
	q2 = (f.V[I_J_Kp2] - f.V[I_J_Kp1])/delta;
	q3 = (f.V[I_J_Kp1] - f.V[I_J_K]  )/delta;
	q4 = (f.V[I_J_K]   - f.V[I_J_Km1])/delta;
	q5 = (f.V[I_J_Km1] - f.V[I_J_Km2])/delta;
}

void ddweno::is()
{
	is1 = tttw*pow(q1-2.0*q2+q3, 2.0) + fourth*pow((q1-4.0*q2+3.0*q3), 2.0);
	is2 = tttw*pow(q2-2.0*q3+q4, 2.0) + fourth*pow((q2-q4), 2.0);
	is3 = tttw*pow(q3-2.0*q4+q5, 2.0) + fourth*pow((3.0*q3-4.0*q4+q5), 2.0);
}

void ddweno::alpha()
{
	alpha1=tenth/pow(epsilon+is1,2.0);
	alpha2=sixten/pow(epsilon+is2,2.0);
	alpha3=treten/pow(epsilon+is3,2.0);
}

void ddweno::weight()
{
	w1=alpha1/(alpha1+alpha2+alpha3);
	w2=alpha2/(alpha1+alpha2+alpha3);
	w3=alpha3/(alpha1+alpha2+alpha3);
}
