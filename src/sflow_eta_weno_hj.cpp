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

#include"sflow_eta_weno_hj.h"
#include"lexer.h"
#include"fdm2D.h"
#include"sflow_flux_HJ_CDS.h"

sflow_eta_weno_hj::sflow_eta_weno_hj(lexer* p):tttw(13.0/12.0),fourth(1.0/4.0),third(1.0/3.0),
			sevsix(7.0/6.0),elvsix(11.0/6.0),sixth(1.0/6.0),fivsix(5.0/6.0),tenth(1.0/10.0),
			sixten(6.0/10.0),treten(3.0/10.0),epsilon(0.000001),smallnum(1.0e-20)
{
        
    pflux = new sflow_flux_HJ_CDS(p);
}

sflow_eta_weno_hj::~sflow_eta_weno_hj()
{
}

void sflow_eta_weno_hj::start(lexer* p, slice& f, int ipol, slice& uvel, slice& vvel, slice &depth,slice &L)
{

    SLICELOOP4
    L(i,j)+=aij(p,f,4,uvel,vvel,depth);

}

double sflow_eta_weno_hj::aij(lexer* p,slice& f,int ipol, slice& uvel, slice& vvel, slice &depth)
{
		pflux->u_flux(4,uvel,iadvec,ivel2);
        pflux->v_flux(4,vvel,jadvec,jvel2);

		L = - iadvec*fx(p,f,ipol,iadvec) 
			- jadvec*fy(p,f,ipol,jadvec);

		return L;
}

double sflow_eta_weno_hj::fx(lexer *p, slice& f, int ipol, double advec)
{
    grad = 0.0;
	
	if(iadvec>0.0)
	{
	iqmin(p,f,ipol);
	is(f);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}
	

	if(iadvec<0.0)
	{
	iqmax(p,f,ipol);
	is(f);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	return grad;
}

double sflow_eta_weno_hj::fy(lexer *p, slice& f, int ipol, double advec)
{
    grad = 0.0;
	
	if(jadvec>0.0)
	{
	jqmin(p,f,ipol);
	is(f);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	
	}
	

	if(jadvec<0.0)
	{
	jqmax(p,f,ipol);
	is(f);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	
	}
	
	return grad;
}


void sflow_eta_weno_hj::iqmin(lexer *p, slice& f, int ipol)
{	
	q1 = (f(i-2,j) - f(i-3,j))/p->DXM;
	q2 = (f(i-1,j) - f(i-2,j))/p->DXM;
	q3 = (f(i,j)   - f(i-1,j))/p->DXM;
	q4 = (f(i+1,j) - f(i,j)  )/p->DXM;
	q5 = (f(i+2,j) - f(i+1,j))/p->DXM;
}

void sflow_eta_weno_hj::jqmin(lexer *p, slice& f, int ipol)
{
	q1 = (f(i,j-2) - f(i,j-3))/p->DXM;
	q2 = (f(i,j-1) - f(i,j-2))/p->DXM;
	q3 = (f(i,j)   - f(i,j-1))/p->DXM;
	q4 = (f(i,j+1) - f(i,j)  )/p->DXM;
	q5 = (f(i,j+2) - f(i,j+1))/p->DXM;
}

void sflow_eta_weno_hj::iqmax(lexer *p, slice& f, int ipol)
{
	q1 = (f(i+3,j) - f(i+2,j))/p->DXM;
	q2 = (f(i+2,j) - f(i+1,j))/p->DXM;
	q3 = (f(i+1,j) - f(i,j)  )/p->DXM;
	q4 = (f(i,j)   - f(i-1,j))/p->DXM;
	q5 = (f(i-1,j) - f(i-2,j))/p->DXM;
}

void sflow_eta_weno_hj::jqmax(lexer *p, slice& f, int ipol)
{
	q1 = (f(i,j+3) - f(i,j+2))/p->DXM;
	q2 = (f(i,j+2) - f(i,j+1))/p->DXM;
	q3 = (f(i,j+1) - f(i,j)  )/p->DXM;
	q4 = (f(i,j)   - f(i,j-1))/p->DXM;
	q5 = (f(i,j-1) - f(i,j-2))/p->DXM;
}

void sflow_eta_weno_hj::is(slice& f)
{
	is1 = tttw*pow(q1 - 2.0*q2 + q3, 2.0) + fourth*pow(q1 - 4.0*q2 + 3.0*q3, 2.0);
	is2 = tttw*pow(q2 - 2.0*q3 + q4, 2.0) + fourth*pow(q2 - q4, 2.0);
	is3 = tttw*pow(q3 - 2.0*q4 + q5, 2.0) + fourth*pow(3.0*q3 - 4.0*q4 + q5, 2.0);
}

void sflow_eta_weno_hj::alpha()
{
	alpha1=tenth/pow(epsilon+is1,2.0);
	alpha2=sixten/pow(epsilon+is2,2.0);
	alpha3=treten/pow(epsilon+is3,2.0);
}

void sflow_eta_weno_hj::weight()
{
	w1=alpha1/(alpha1+alpha2+alpha3);
	w2=alpha2/(alpha1+alpha2+alpha3);
	w3=alpha3/(alpha1+alpha2+alpha3);
}
