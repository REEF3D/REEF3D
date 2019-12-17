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

#include"weno_hj_N.h"
#include"lexer.h"
#include"fdm.h"
#include"flux_HJ_CDS2.h"
#include"flux_HJ_CDS4.h"
#include"flux_HJ_CDS2_vrans.h"
 
weno_hj_N::weno_hj_N(lexer* p):tttw(13.0/12.0),fourth(1.0/4.0),third(1.0/3.0),
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

weno_hj_N::~weno_hj_N()
{
}

void weno_hj_N::start(lexer* p, fdm* a, fieldint& cv, vec &x, int ipol, field& uvel, field& vvel, field& wvel)
{

        sizeM = p->sizeM4;
        

		n=0;
		LOOP
		{
		x.V[n]+=aij(p,a,x,cv,4,uvel,vvel,wvel,a->C4);
		++n;
		}
	
}

double weno_hj_N::aij(lexer* p,fdm* a,vec& x,fieldint& cv,int ipol, field& uvel, field& vvel, field& wvel, cpt &C)
{

		pflux->u_flux(a,ipol,uvel,iadvec,ivel2);
        pflux->v_flux(a,ipol,vvel,jadvec,jvel2);
        pflux->w_flux(a,ipol,wvel,kadvec,kvel2);
		
		L = -iadvec*ddx(p,a,x,C) - jadvec*ddy(p,a,x,C) - kadvec*ddz(p,a,x,C);

		return L;
}

double weno_hj_N::ddx(lexer* p,fdm* a, vec& x, cpt &C)
{
    grad = 0.0;
	
	if(iadvec>0.0)
	{
	iqmin(a,x,p->dx,C);
	is();
	alpha();
	weight();
    
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}
	

	if(iadvec<0.0)
	{
	iqmax(a,x,p->dx,C);
	is();
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	return grad;
}

double weno_hj_N::ddy(lexer* p,fdm* a, vec& x, cpt &C)
{
    grad = 0.0;
	
	if(jadvec>0.0)
	{
	jqmin(a,x,p->dx,C);
	is();
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	
	}
	

	if(jadvec<0.0)
	{
	jqmax(a,x,p->dx,C);
	is();
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	
	}
	
	return grad;
}

double weno_hj_N::ddz(lexer* p,fdm* a, vec& x, cpt &C)
{
    grad = 0.0;
	
	if(kadvec>0.0)
	{
	kqmin(a,x,p->dx,C);
	is();
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}
	 
	
	if(kadvec<0.0)
	{
	kqmax(a,x,p->dx,C);
	is();
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}
	
	return grad;
}

void weno_hj_N::iqmin(fdm *a, vec& x, double dx, cpt &C)
{	
	q1 = (x.V[Im2_J_K] - x.V[Im3_J_K])/dx;
	q2 = (x.V[Im1_J_K] - x.V[Im2_J_K])/dx;
	q3 = (x.V[I_J_K]   - x.V[Im1_J_K])/dx;
	q4 = (x.V[Ip1_J_K] - x.V[I_J_K]  )/dx;
	q5 = (x.V[Ip2_J_K] - x.V[Ip1_J_K])/dx;
}

void weno_hj_N::jqmin(fdm *a, vec& x, double dx, cpt &C)
{	
	q1 = (x.V[I_Jm2_K] - x.V[I_Jm3_K])/dx;
	q2 = (x.V[I_Jm1_K] - x.V[I_Jm2_K])/dx;
	q3 = (x.V[I_J_K]   - x.V[I_Jm1_K])/dx;
	q4 = (x.V[I_Jp1_K] - x.V[I_J_K]  )/dx;
	q5 = (x.V[I_Jp2_K] - x.V[I_Jp1_K])/dx;
}

void weno_hj_N::kqmin(fdm *a, vec& x, double dx, cpt &C)
{	
	q1 = (x.V[I_J_Km2] - x.V[I_J_Km3])/dx;
	q2 = (x.V[I_J_Km1] - x.V[I_J_Km2])/dx;
	q3 = (x.V[I_J_K]   - x.V[I_J_Km1])/dx;
	q4 = (x.V[I_J_Kp1] - x.V[I_J_K]  )/dx;
	q5 = (x.V[I_J_Kp2] - x.V[I_J_Kp1])/dx;
}

void weno_hj_N::iqmax(fdm *a, vec& x, double dx, cpt &C)
{	
	q1 = (x.V[Ip3_J_K] - x.V[Ip2_J_K])/dx;
	q2 = (x.V[Ip2_J_K] - x.V[Ip1_J_K])/dx;
	q3 = (x.V[Ip1_J_K] - x.V[I_J_K]  )/dx;
	q4 = (x.V[I_J_K]   - x.V[Im1_J_K])/dx;
	q5 = (x.V[Im1_J_K] - x.V[Im2_J_K])/dx;
}

void weno_hj_N::jqmax(fdm *a, vec& x, double dx, cpt &C)
{	
	q1 = (x.V[I_Jp3_K] - x.V[I_Jp2_K])/dx;
	q2 = (x.V[I_Jp2_K] - x.V[I_Jp1_K])/dx;
	q3 = (x.V[I_Jp1_K] - x.V[I_J_K]  )/dx;
	q4 = (x.V[I_J_K]   - x.V[I_Jm1_K])/dx;
	q5 = (x.V[I_Jm1_K] - x.V[I_Jm2_K])/dx;
}

void weno_hj_N::kqmax(fdm *a, vec& x, double dx, cpt &C)
{
	q1 = (x.V[I_J_Kp3] - x.V[I_J_Kp2])/dx;
	q2 = (x.V[I_J_Kp2] - x.V[I_J_Kp1])/dx;
	q3 = (x.V[I_J_Kp1] - x.V[I_J_K]  )/dx;
	q4 = (x.V[I_J_K]   - x.V[I_J_Km1])/dx;
	q5 = (x.V[I_J_Km1] - x.V[I_J_Km2])/dx;
}

void weno_hj_N::is()
{
	is1 = tttw*pow(q1-2.0*q2+q3, 2.0) + fourth*pow(q1-4.0*q2+3.0*q3, 2.0);
	is2 = tttw*pow(q2-2.0*q3+q4, 2.0) + fourth*pow(q2-q4, 2.0);
	is3 = tttw*pow(q3-2.0*q4+q5, 2.0) + fourth*pow(3.0*q3-4.0*q4+q5, 2.0);
}

void weno_hj_N::alpha()
{
	alpha1=tenth/pow(epsilon+is1,2.0);
	alpha2=sixten/pow(epsilon+is2,2.0);
	alpha3=treten/pow(epsilon+is3,2.0);
}

void weno_hj_N::weight()
{
	w1=alpha1/(alpha1+alpha2+alpha3);
	w2=alpha2/(alpha1+alpha2+alpha3);
	w3=alpha3/(alpha1+alpha2+alpha3);
}

