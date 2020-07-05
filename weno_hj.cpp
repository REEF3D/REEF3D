/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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
--------------------------------------------------------------------*/

#include"weno_hj.h"
#include"lexer.h"
#include"fdm.h"
#include"flux_HJ_CDS2.h"
#include"flux_HJ_CDS4.h"
#include"flux_HJ_CDS2_vrans.h"

weno_hj::weno_hj(lexer* p):tttw(13.0/12.0),fourth(1.0/4.0),third(1.0/3.0),
			sevsix(7.0/6.0),elvsix(11.0/6.0),sixth(1.0/6.0),fivsix(5.0/6.0),tenth(1.0/10.0),
			sixten(6.0/10.0),treten(3.0/10.0),epsilon(0.000001),smallnum(1.0e-20)
{
    if(p->B269==0 && p->D11!=4)
    pflux = new flux_HJ_CDS2(p);
    
    if(p->B269==0 && p->D11==4)
    pflux = new flux_HJ_CDS4(p);
    
    if(p->B269>=1 || p->S10==2)
    pflux = new flux_HJ_CDS2_vrans(p);
}

weno_hj::~weno_hj()
{
}

void weno_hj::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{
	
	
    if(ipol==1)
	{
		n=0;
		ULOOP
		{
		a->F(i,j,k)+=aij(p,a,b,1,uvel,vvel,wvel);
		++n;
		}
	}

    if(ipol==2)
	{
		n=0;
		VLOOP
		{
		a->G(i,j,k)+=aij(p,a,b,2,uvel,vvel,wvel);
		++n;
		}
	}

    if(ipol==3)
	{
		n=0;
		WLOOP
		{
		a->H(i,j,k)+=aij(p,a,b,3,uvel,vvel,wvel);
		++n;
		}
	}

    if(ipol==4)
	{
		n=0;
		FLUIDLOOP
		{
		a->L(i,j,k)+=aij(p,a,b,4,uvel,vvel,wvel);
		++n;
		}
	}
}

double weno_hj::aij(lexer* p,fdm* a,field& b,int ipol, field& uvel, field& vvel, field& wvel)
{

		pflux->u_flux(a,ipol,uvel,iadvec,ivel2);
        pflux->v_flux(a,ipol,vvel,jadvec,jvel2);
        pflux->w_flux(a,ipol,wvel,kadvec,kvel2);
		
		L = -iadvec*ddx(p,a,b) - jadvec*ddy(p,a,b) - kadvec*ddz(p,a,b);

		return L;
}

double weno_hj::ddx(lexer* p,fdm* a, field& b)
{
    grad = 0.0;
	
	if(iadvec>0.0)
	{
	iqmin(a,b,p->DXM);
	is();
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}
	

	if(iadvec<0.0)
	{
	iqmax(a,b,p->DXM);
	is();
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	return grad;
}

double weno_hj::ddy(lexer* p,fdm* a, field& b)
{
    grad = 0.0;
	
	if(jadvec>0.0)
	{
	jqmin(a,b,p->DXM);
	is();
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	
	}
	

	if(jadvec<0.0)
	{
	jqmax(a,b,p->DXM);
	is();
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	
	}
	
	return grad;
}

double weno_hj::ddz(lexer* p,fdm* a, field& b)
{
    grad = 0.0;
	
	if(kadvec>0.0)
	{
	kqmin(a,b,p->DXM);
	is();
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}
	 
	
	if(kadvec<0.0)
	{
	kqmax(a,b,p->DXM);
	is();
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}
	
	return grad;
}

void weno_hj::iqmin(fdm *a, field& f, double dx)
{	
	q1 = (f(i-2,j,k) - f(i-3,j,k))/dx;
	q2 = (f(i-1,j,k) - f(i-2,j,k))/dx;
	q3 = (f(i,j,k)   - f(i-1,j,k))/dx;
	q4 = (f(i+1,j,k) - f(i,j,k)  )/dx;
	q5 = (f(i+2,j,k) - f(i+1,j,k))/dx;
}

void weno_hj::jqmin(fdm *a, field& f, double dx)
{	
	q1 = (f(i,j-2,k) - f(i,j-3,k))/dx;
	q2 = (f(i,j-1,k) - f(i,j-2,k))/dx;
	q3 = (f(i,j,k)   - f(i,j-1,k))/dx;
	q4 = (f(i,j+1,k) - f(i,j,k)  )/dx;
	q5 = (f(i,j+2,k) - f(i,j+1,k))/dx;
}

void weno_hj::kqmin(fdm *a, field& f, double dx)
{	
	q1 = (f(i,j,k-2) - f(i,j,k-3))/dx;
	q2 = (f(i,j,k-1) - f(i,j,k-2))/dx;
	q3 = (f(i,j,k)   - f(i,j,k-1))/dx;
	q4 = (f(i,j,k+1) - f(i,j,k)  )/dx;
	q5 = (f(i,j,k+2) - f(i,j,k+1))/dx;
}

void weno_hj::iqmax(fdm *a, field& f, double dx)
{	
	q1 = (f(i+3,j,k) - f(i+2,j,k))/dx;
	q2 = (f(i+2,j,k) - f(i+1,j,k))/dx;
	q3 = (f(i+1,j,k) - f(i,j,k)  )/dx;
	q4 = (f(i,j,k)   - f(i-1,j,k))/dx;
	q5 = (f(i-1,j,k) - f(i-2,j,k))/dx;
}

void weno_hj::jqmax(fdm *a, field& f, double dx)
{	
	q1 = (f(i,j+3,k) - f(i,j+2,k))/dx;
	q2 = (f(i,j+2,k) - f(i,j+1,k))/dx;
	q3 = (f(i,j+1,k) - f(i,j,k)  )/dx;
	q4 = (f(i,j,k)   - f(i,j-1,k))/dx;
	q5 = (f(i,j-1,k) - f(i,j-2,k))/dx;
}

void weno_hj::kqmax(fdm *a, field& f, double dx)
{
	q1 = (f(i,j,k+3) - f(i,j,k+2))/dx;
	q2 = (f(i,j,k+2) - f(i,j,k+1))/dx;
	q3 = (f(i,j,k+1) - f(i,j,k)  )/dx;
	q4 = (f(i,j,k)   - f(i,j,k-1))/dx;
	q5 = (f(i,j,k-1) - f(i,j,k-2))/dx;
}

void weno_hj::is()
{
	is1 = tttw*pow(q1-2.0*q2+q3, 2.0) + fourth*pow(q1-4.0*q2+3.0*q3, 2.0);
	is2 = tttw*pow(q2-2.0*q3+q4, 2.0) + fourth*pow(q2-q4, 2.0);
	is3 = tttw*pow(q3-2.0*q4+q5, 2.0) + fourth*pow(3.0*q3-4.0*q4+q5, 2.0);
}

void weno_hj::alpha()
{
	alpha1=tenth/pow(epsilon+is1,2.0);
	alpha2=sixten/pow(epsilon+is2,2.0);
	alpha3=treten/pow(epsilon+is3,2.0);
}

void weno_hj::weight()
{
	w1=alpha1/(alpha1+alpha2+alpha3);
	w2=alpha2/(alpha1+alpha2+alpha3);
	w3=alpha3/(alpha1+alpha2+alpha3);
}

