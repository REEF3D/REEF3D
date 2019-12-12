/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"weno_flux_N.h"
#include"lexer.h"
#include"fdm.h"
#include"flux_face_CDS2.h"
#include"flux_face_CDS4.h"
#include"flux_face_CDS2_vrans.h"
#include"flux_face_FOU.h"
#include"flux_face_FOU_vrans.h"
#include"flux_face_QOU.h"

weno_flux_N::weno_flux_N(lexer* p):tttw(13.0/12.0),fourth(1.0/4.0),third(1.0/3.0),
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

weno_flux_N::~weno_flux_N()
{
}

void weno_flux_N::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{   
    
    
    if(ipol==1)
	{
	fillxvec1(p,a,b);
        n=0;
        ULOOP
        {
        a->F(i,j,k)+=aij(p,a,b,1,uvel,vvel,wvel,a->C1);
        ++n;
        }
	}

    if(ipol==2)
	{
	fillxvec2(p,a,b);
        n=0;
        VLOOP
        {
        a->G(i,j,k)+=aij(p,a,b,2,uvel,vvel,wvel,a->C2);
        ++n;
        }
	}

    if(ipol==3)
	{
	fillxvec3(p,a,b);
        n=0;
        WLOOP
        {
        a->H(i,j,k)+=aij(p,a,b,3,uvel,vvel,wvel,a->C3);
        ++n;
        }
	}

    if(ipol==4)
	{
	fillxvec4(p,a,b);
        n=0;
        LOOP
        {
        a->L(i,j,k)+=aij(p,a,b,4,uvel,vvel,wvel,a->C4);
        ++n;
        }
	}
}

double weno_flux_N::aij(lexer* p,fdm* a,field& b,int ipol, field& uvel, field& vvel, field& wvel, cpt &C)
{
        n0=n;
        
		pflux->u_flux(a,ipol,uvel,ivel1,ivel2);
        pflux->v_flux(a,ipol,vvel,jvel1,jvel2);
        pflux->w_flux(a,ipol,wvel,kvel1,kvel2);
		
		n=Im1_J_K;
		fu1 = fx(p,a,b,uvel,ivel1,C);
        n=n0;
		
		fu2 = fx(p,a,b,uvel,ivel2,C);


		n=I_Jm1_K;
		fv1 = fy(p,a,b,vvel,jvel1,C);
        n=n0;
		
		fv2 = fy(p,a,b,vvel,jvel2,C);


		n=I_J_Km1;
		fw1 = fz(p,a,b,wvel,kvel1,C);
        n=n0;
		
		fw2 = fz(p,a,b,wvel,kvel2,C);
		
		
		L =   - ((ivel2*fu2-ivel1*fu1)/p->dx) 
		      - ((jvel2*fv2-jvel1*fv1)/p->dx) 
			  - ((kvel2*fw2-kvel1*fw1)/p->dx);
			  
			  
		return L;
}

double weno_flux_N::fx(lexer *p,fdm *a, field& b, field& uvel, double advec, cpt &C)
{
    grad = 0.0;

	if(advec>0.0)
	{
	iqmin(p,a,b,uvel,C);
	is(b);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	if(advec<0.0)
	{
	iqmax(p,a,b,uvel,C);
	is(b);
	alpha();
	weight();
   
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	return grad;
}

double weno_flux_N::fy(lexer *p,fdm *a, field& b, field& vvel, double advec, cpt &C)
{
    grad = 0.0;

	if(advec>0.0)
	{
	jqmin(p,a,b,vvel,C);
	is(b);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	if(advec<0.0)
	{
	jqmax(p,a,b,vvel,C);
	is(b);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}
	
	return grad;
}

double weno_flux_N::fz(lexer *p,fdm *a, field& b, field& wvel, double advec, cpt &C)
{
    grad = 0.0;

	if(advec>0.0)
	{
	kqmin(p,a,b,wvel,C);
	is(b);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	if(advec<0.0)
	{
	kqmax(p,a,b,wvel,C);
	is(b);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	return grad;
}

void weno_flux_N::iqmin(lexer *p,fdm *a, field& f, field& uvel, cpt &C)
{		
	q1 = a->xvec.V[Im2_J_K];
	q2 = a->xvec.V[Im1_J_K];
	q3 = a->xvec.V[I_J_K];
	q4 = a->xvec.V[Ip1_J_K];
	q5 = a->xvec.V[Ip2_J_K];
}

void weno_flux_N::jqmin(lexer *p,fdm *a, field& f, field& vvel, cpt &C)
{
	q1 = a->xvec.V[I_Jm2_K];
	q2 = a->xvec.V[I_Jm1_K];
	q3 = a->xvec.V[I_J_K];
	q4 = a->xvec.V[I_Jp1_K];
	q5 = a->xvec.V[I_Jp2_K];
}

void weno_flux_N::kqmin(lexer *p,fdm *a, field& f, field& wvel, cpt &C)
{
	q1 = a->xvec.V[I_J_Km2];
	q2 = a->xvec.V[I_J_Km1];
	q3 = a->xvec.V[I_J_K];
	q4 = a->xvec.V[I_J_Kp1];
	q5 = a->xvec.V[I_J_Kp2];
}

void weno_flux_N::iqmax(lexer *p,fdm *a, field& f, field& uvel, cpt &C)
{	
	q1 = a->xvec.V[Ip3_J_K];
	q2 = a->xvec.V[Ip2_J_K];
	q3 = a->xvec.V[Ip1_J_K];
	q4 = a->xvec.V[I_J_K];
	q5 = a->xvec.V[Im1_J_K];
}

void weno_flux_N::jqmax(lexer *p,fdm *a, field& f, field& vvel, cpt &C)
{
	q1 = a->xvec.V[I_Jp3_K];
	q2 = a->xvec.V[I_Jp2_K];
	q3 = a->xvec.V[I_Jp1_K];
	q4 = a->xvec.V[I_J_K];
	q5 = a->xvec.V[I_Jm1_K];
}

void weno_flux_N::kqmax(lexer *p,fdm *a, field& f, field& wvel, cpt &C)
{
	q1 = a->xvec.V[I_J_Kp3];
	q2 = a->xvec.V[I_J_Kp2];
	q3 = a->xvec.V[I_J_Kp1];
	q4 = a->xvec.V[I_J_K];
	q5 = a->xvec.V[I_J_Km1];
}

void weno_flux_N::is(field& f)
{
	is1 = tttw*pow(q1 - 2.0*q2 + q3, 2.0) + fourth*pow(q1 - 4.0*q2 + 3.0*q3, 2.0);
	is2 = tttw*pow(q2 - 2.0*q3 + q4, 2.0) + fourth*pow(q2 - q4, 2.0);
	is3 = tttw*pow(q3 - 2.0*q4 + q5, 2.0) + fourth*pow(3.0*q3 - 4.0*q4 + q5, 2.0);
}

void weno_flux_N::alpha()
{
	alpha1=tenth/pow(epsilon+is1,2.0);
	alpha2=sixten/pow(epsilon+is2,2.0);
	alpha3=treten/pow(epsilon+is3,2.0);
}

void weno_flux_N::weight()
{
	w1=alpha1/(alpha1+alpha2+alpha3);
	w2=alpha2/(alpha1+alpha2+alpha3);
	w3=alpha3/(alpha1+alpha2+alpha3);
}
