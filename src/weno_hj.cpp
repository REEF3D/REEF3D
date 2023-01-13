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
		a->F(i,j,k)+=aij(p,a,b,1,uvel,vvel,wvel,p->DRDXN,p->DSDYP,p->DTDZP);
		++n;
		}
	}

    if(ipol==2)
	{
		n=0;
		VLOOP
		{
		a->G(i,j,k)+=aij(p,a,b,2,uvel,vvel,wvel,p->DRDXP,p->DSDYN,p->DTDZP);
		++n;
		}
	}

    if(ipol==3)
	{
		n=0;
		WLOOP
		{
		a->H(i,j,k)+=aij(p,a,b,3,uvel,vvel,wvel,p->DRDXP,p->DSDYP,p->DTDZN);
		++n;
		}
	}

    if(ipol==4)
	{
		n=0;
		FLUIDLOOP
		{
		a->L(i,j,k)+=aij(p,a,b,4,uvel,vvel,wvel,p->DRDXP,p->DSDYP,p->DTDZP);
		++n;
		}
	}
}

double weno_hj::aij(lexer* p,fdm* a,field& b,int ipol, field& uvel, field& vvel, field& wvel, double *DRDX, double *DSDY, double *DTDZ)
{

		pflux->u_flux(a,ipol,uvel,iadvec,ivel2);
        pflux->v_flux(a,ipol,vvel,jadvec,jvel2);
        pflux->w_flux(a,ipol,wvel,kadvec,kvel2);
		
		L = -iadvec*ddx(p,a,b)*DRDX[IP]  - jadvec*ddy(p,a,b)*DSDY[JP]  - kadvec*ddz(p,a,b)*DTDZ[KP] ;

		return L;
}

double weno_hj::ddx(lexer* p,fdm* a, field& b)
{
    grad = 0.0;
	
	if(iadvec>0.0)
	{
	iqmin(b,p->DRM,p->DRDXN);
	is();
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}
	

	if(iadvec<0.0)
	{
	iqmax(b,p->DRM,p->DRDXN);
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
	jqmin(b,p->DSM,p->DSDYN);
	is();
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	
	}
	

	if(jadvec<0.0)
	{
	jqmax(b,p->DSM,p->DSDYN);
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
	kqmin(b,p->DTM,p->DTDZN);
	is();
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}
	 
	
	if(kadvec<0.0)
	{
	kqmax(b,p->DTM,p->DTDZN);
	is();
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}
	
	return grad;
}

void weno_hj::iqmin(field& f, double DRM, double *DRDX)
{	
	q1 = (f(i-2,j,k) - f(i-3,j,k))/DRM;
	q2 = (f(i-1,j,k) - f(i-2,j,k))/DRM;
	q3 = (f(i,j,k)   - f(i-1,j,k))/DRM;
	q4 = (f(i+1,j,k) - f(i,j,k)  )/DRM;
	q5 = (f(i+2,j,k) - f(i+1,j,k))/DRM;
}

void weno_hj::jqmin(field& f, double DSM, double *DSDY)
{	
	q1 = (f(i,j-2,k) - f(i,j-3,k))/DSM;
	q2 = (f(i,j-1,k) - f(i,j-2,k))/DSM;
	q3 = (f(i,j,k)   - f(i,j-1,k))/DSM;
	q4 = (f(i,j+1,k) - f(i,j,k)  )/DSM;
	q5 = (f(i,j+2,k) - f(i,j+1,k))/DSM;
}

void weno_hj::kqmin(field& f, double DTM, double *DTDZ)
{	
	q1 = (f(i,j,k-2) - f(i,j,k-3))/DTM;
	q2 = (f(i,j,k-1) - f(i,j,k-2))/DTM;
	q3 = (f(i,j,k)   - f(i,j,k-1))/DTM;
	q4 = (f(i,j,k+1) - f(i,j,k)  )/DTM;
	q5 = (f(i,j,k+2) - f(i,j,k+1))/DTM;
}

void weno_hj::iqmax(field& f, double DRM, double *DRDX)
{	
	q1 = (f(i+3,j,k) - f(i+2,j,k))/DRM;
	q2 = (f(i+2,j,k) - f(i+1,j,k))/DRM;
	q3 = (f(i+1,j,k) - f(i,j,k)  )/DRM;
	q4 = (f(i,j,k)   - f(i-1,j,k))/DRM;
	q5 = (f(i-1,j,k) - f(i-2,j,k))/DRM;
}

void weno_hj::jqmax(field& f, double DSM, double *DSDY)
{	
	q1 = (f(i,j+3,k) - f(i,j+2,k))/DSM;
	q2 = (f(i,j+2,k) - f(i,j+1,k))/DSM;
	q3 = (f(i,j+1,k) - f(i,j,k)  )/DSM;
	q4 = (f(i,j,k)   - f(i,j-1,k))/DSM;
	q5 = (f(i,j-1,k) - f(i,j-2,k))/DSM;
}

void weno_hj::kqmax(field& f, double DTM, double *DTDZ)
{
	q1 = (f(i,j,k+3) - f(i,j,k+2))/DTM;
	q2 = (f(i,j,k+2) - f(i,j,k+1))/DTM;
	q3 = (f(i,j,k+1) - f(i,j,k)  )/DTM;
	q4 = (f(i,j,k)   - f(i,j,k-1))/DTM;
	q5 = (f(i,j,k-1) - f(i,j,k-2))/DTM;
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

