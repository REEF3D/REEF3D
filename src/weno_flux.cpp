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

#include"weno_flux.h"
#include"lexer.h"
#include"fdm.h"
#include"flux_face_CDS2.h"
#include"flux_face_CDS4.h"
#include"flux_face_CDS2_vrans.h"
#include"flux_face_FOU.h"
#include"flux_face_FOU_vrans.h"
#include"flux_face_QOU.h"

weno_flux::weno_flux(lexer* p):tttw(13.0/12.0),fourth(1.0/4.0),third(1.0/3.0),
			sevsix(7.0/6.0),elvsix(11.0/6.0),sixth(1.0/6.0),fivsix(5.0/6.0),tenth(1.0/10.0),
			sixten(6.0/10.0),treten(3.0/10.0),epsilon(0.000001),smallnum(1.0e-20)
{
    
    if(p->B269==0)
    {
        if(p->D11==1)
        pflux = new flux_face_FOU(p);
        
        if(p->D11==2)
        pflux = new flux_face_CDS2(p);
        
        if(p->D11==3)
        pflux = new flux_face_QOU(p);
        
        if(p->D11==4)
        pflux = new flux_face_CDS4(p);
    }
    
    if(p->B269>=1 || p->S10==2)
    {
        if(p->D11==1)
        pflux = new flux_face_FOU_vrans(p);
        
        if(p->D11==2)
        pflux = new flux_face_CDS2_vrans(p);
        
        if(p->D11==3)
        pflux = new flux_face_FOU_vrans(p);
        
        if(p->D11==4)
        pflux = new flux_face_CDS2(p);
    }
}

weno_flux::~weno_flux()
{
}

void weno_flux::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{
    if(ipol==1)
    ULOOP
    a->F(i,j,k)+=aij(p,a,b,1,uvel,vvel,wvel,p->DRDXN,p->DSDYP,p->DTDZP);

    if(ipol==2)
    VLOOP
    a->G(i,j,k)+=aij(p,a,b,2,uvel,vvel,wvel,p->DRDXP,p->DSDYN,p->DTDZP);

    if(ipol==3)
    WLOOP
    a->H(i,j,k)+=aij(p,a,b,3,uvel,vvel,wvel,p->DRDXP,p->DSDYP,p->DTDZN);

    if(ipol==4)
    LOOP
    a->L(i,j,k)+=aij(p,a,b,4,uvel,vvel,wvel,p->DRDXP,p->DSDYP,p->DTDZP);

    if(ipol==5)
    LOOP
    a->L(i,j,k)+=aij(p,a,b,5,uvel,vvel,wvel,p->DRDXP,p->DSDYP,p->DTDZP);

}

double weno_flux::aij(lexer* p,fdm* a,field& b,int ipol, field& uvel, field& vvel, field& wvel, double *DRDX, double *DSDY, double *DTDZ)
{

		pflux->u_flux(a,ipol,uvel,ivel1,ivel2);
        pflux->v_flux(a,ipol,vvel,jvel1,jvel2);
        pflux->w_flux(a,ipol,wvel,kvel1,kvel2);

		
		i-=1;
		fu1 = fx(p,a,b,uvel,ipol,ivel1);
		i+=1;
		
		fu2 = fx(p,a,b,uvel,ipol,ivel2);


		
		j-=1;
		fv1 = fy(p,a,b,vvel,ipol,jvel1);
		j+=1;
		
		fv2 = fy(p,a,b,vvel,ipol,jvel2);



		k-=1;
		fw1 = fz(p,a,b,wvel,ipol,kvel1);
		k+=1;
		
		fw2 = fz(p,a,b,wvel,ipol,kvel2);
		
		
		L =   - ((ivel2*fu2-ivel1*fu1)/p->DRM)*DRDX[IP] 
		      - ((jvel2*fv2-jvel1*fv1)/p->DSM)*DSDY[JP] 
			  - ((kvel2*fw2-kvel1*fw1)/p->DTM)*DTDZ[KP];
			  
			  
		return L;
}

double weno_flux::fx(lexer *p,fdm *a, field& b, field& uvel, int ipol, double advec)
{
    grad = 0.0;

	if(advec>0.0)
	{
	iqmin(p,a,b,uvel,ipol);
	is(b);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
          
	}

	if(advec<0.0)
	{
	iqmax(p,a,b,uvel,ipol);
	is(b);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	return grad;
}

double weno_flux::fy(lexer *p,fdm *a, field& b, field& vvel, int ipol, double advec)
{
    grad = 0.0;

	if(advec>0.0)
	{
	jqmin(p,a,b,vvel,ipol);
	is(b);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	if(advec<0.0)
	{
	jqmax(p,a,b,vvel,ipol);
	is(b);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}
	
	return grad;
}

double weno_flux::fz(lexer *p,fdm *a, field& b, field& wvel, int ipol, double advec)
{
    grad = 0.0;

	if(advec>0.0)
	{
	kqmin(p,a,b,wvel,ipol);
	is(b);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	if(advec<0.0)
	{
	kqmax(p,a,b,wvel,ipol);
	is(b);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	return grad;
}

void weno_flux::iqmin(lexer *p,fdm *a, field& f, field& uvel, int ipol)
{	
	q1 = f(i-2,j,k);
	q2 = f(i-1,j,k);
	q3 = f(i,j,k);
	q4 = f(i+1,j,k);
	q5 = f(i+2,j,k);
}

void weno_flux::jqmin(lexer *p,fdm *a, field& f, field& vvel, int ipol)
{
	q1 = f(i,j-2,k);
	q2 = f(i,j-1,k);
	q3 = f(i,j,k);
	q4 = f(i,j+1,k);
	q5 = f(i,j+2,k);
}

void weno_flux::kqmin(lexer *p,fdm *a, field& f, field& wvel, int ipol)
{
	q1 = f(i,j,k-2);
	q2 = f(i,j,k-1);
	q3 = f(i,j,k);
	q4 = f(i,j,k+1);
	q5 = f(i,j,k+2);
}

void weno_flux::iqmax(lexer *p,fdm *a, field& f, field& uvel, int ipol)
{
	q1 = f(i+3,j,k);
	q2 = f(i+2,j,k);
	q3 = f(i+1,j,k);
	q4 = f(i,j,k);
	q5 = f(i-1,j,k);
}

void weno_flux::jqmax(lexer *p,fdm *a, field& f, field& vvel, int ipol)
{
	q1 = f(i,j+3,k);
	q2 = f(i,j+2,k);
	q3 = f(i,j+1,k);
	q4 = f(i,j,k);
	q5 = f(i,j-1,k);
}

void weno_flux::kqmax(lexer *p,fdm *a, field& f, field& wvel, int ipol)
{
	q1 = f(i,j,k+3);
	q2 = f(i,j,k+2);
	q3 = f(i,j,k+1);
	q4 = f(i,j,k);
	q5 = f(i,j,k-1);
}

void weno_flux::is(field& f)
{
	is1 = tttw*pow(q1 - 2.0*q2 + q3, 2.0) + fourth*pow(q1 - 4.0*q2 + 3.0*q3, 2.0);
	is2 = tttw*pow(q2 - 2.0*q3 + q4, 2.0) + fourth*pow(q2 - q4, 2.0);
	is3 = tttw*pow(q3 - 2.0*q4 + q5, 2.0) + fourth*pow(3.0*q3 - 4.0*q4 + q5, 2.0);
}

void weno_flux::alpha()
{
	alpha1=tenth/pow(epsilon+is1,2.0);
	alpha2=sixten/pow(epsilon+is2,2.0);
	alpha3=treten/pow(epsilon+is3,2.0);
}

void weno_flux::weight()
{
	w1=alpha1/(alpha1+alpha2+alpha3);
	w2=alpha2/(alpha1+alpha2+alpha3);
	w3=alpha3/(alpha1+alpha2+alpha3);
}
