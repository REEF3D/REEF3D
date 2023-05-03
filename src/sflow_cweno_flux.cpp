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
--------------------------------------------------------------------*/

#include"sflow_cweno_flux.h"
#include"lexer.h"
#include"fdm2D.h"
#include"sflow_flux_face_FOU.h"
#include"sflow_flux_face_CDS.h"
#include"sflow_flux_HJ_CDS.h"
#include"sflow_flux_face_C_FOU.h"
#include"sflow_flux_face_C_CDS.h"
#include"sflow_flux_face_C_HJ.h"

#define HXIJ (fabs(b->hx(i,j))>1.0e-20?b->hx(i,j):1.0e20)
#define HYIJ (fabs(b->hy(i,j))>1.0e-20?b->hy(i,j):1.0e20)
#define HP (fabs(b->hp(i,j))>1.0e-20?b->hp(i,j):1.0e20)

sflow_cweno_flux::sflow_cweno_flux(lexer* p, fdm2D *b):tttw(13.0/12.0),fourth(1.0/4.0),third(1.0/3.0),
			sevsix(7.0/6.0),elvsix(11.0/6.0),sixth(1.0/6.0),fivsix(5.0/6.0),tenth(1.0/10.0),
			sixten(6.0/10.0),treten(3.0/10.0),epsilon(0.000001),smallnum(1.0e-20)
{
    if(p->A216==1)
    pflux = new sflow_flux_face_C_FOU(p,b);
        
    if(p->A216==2)
    pflux = new sflow_flux_face_C_CDS(p,b);
        
    if(p->A216==4)
    pflux = new sflow_flux_face_C_HJ(p,b);
}

sflow_cweno_flux::~sflow_cweno_flux()
{
}

void sflow_cweno_flux::start(lexer* p, fdm2D* b, slice& f, int ipol, slice& uvel, slice& vvel)
{
    if(ipol==1)
    SLICELOOP1
    b->F(i,j)+=aij(p,b,f,1,uvel,vvel);

    if(ipol==2)
    SLICELOOP2
    b->G(i,j)+=aij(p,b,f,2,uvel,vvel);
    
    if(ipol==4)
    SLICELOOP4
    b->L(i,j)+=aij(p,b,f,4,uvel,vvel);

    if(ipol==5)
    SLICELOOP4
    b->L(i,j)+=aij(p,b,f,5,uvel,vvel);

}

double sflow_cweno_flux::aij(lexer* p,fdm2D* b,slice& f,int ipol, slice& uvel, slice& vvel)
{
		pflux->u_flux(ipol,uvel,ivel1,ivel2);
        pflux->v_flux(ipol,vvel,jvel1,jvel2);
        

		i-=1;
		fu1 = fx(p,b,f,ipol,ivel1);
		i+=1;
		
		fu2 = fx(p,b,f,ipol,ivel2);

		
		j-=1;
		fv1 = fy(p,b,f,ipol,jvel1);
		j+=1;
		
		fv2 = fy(p,b,f,ipol,jvel2);
        
        
		
		// X
		dx = (ivel2*fu2-ivel1*fu1)/p->DXM; 
        
            if(ipol<=2)
            dx -= f(i,j)*(ivel2-ivel1)/p->DXM;
        
            if(ipol==1)
            dx/=HXIJ;
            
            if(ipol==2)
            dx/=HYIJ;
            
            if(ipol==4)
            dx/=HP;
        
        // Y
		dy = (jvel2*fv2-jvel1*fv1)/p->DXM;
        
            if(ipol<=2)
            dy -= f(i,j)*(jvel2-jvel1)/p->DXM;
            
            if(ipol==1)
            dy/=HXIJ;
            
            if(ipol==2)
            dy/=HYIJ;
        
            if(ipol==4)
            dy/=HP;
        
        L = -dx-dy;
  			  
		return L;
}

double sflow_cweno_flux::fx(lexer *p,fdm2D *b, slice& f, int ipol, double advec)
{
    grad = 0.0;

	if(advec>0.0)
	{
	iqmin(p,b,f,ipol);
	is_min(b,f);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	if(advec<0.0)
	{
	iqmax(p,b,f,ipol);
	is_max(b,f);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	return grad;
}

double sflow_cweno_flux::fy(lexer *p,fdm2D *b, slice& f, int ipol, double advec)
{
    grad = 0.0;

	if(advec>0.0)
	{
	jqmin(p,b,f,ipol);
	is_min(b,f);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	if(advec<0.0)
	{
	jqmax(p,b,f,ipol);
	is_max(b,f);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}
	
	return grad;
}

void sflow_cweno_flux::iqmin(lexer *p,fdm2D *b, slice& f, int ipol)
{	
	q1 = f(i-2,j);
	q2 = f(i-1,j);
	q3 = f(i,j);
	q4 = f(i+1,j);
	q5 = f(i+2,j);
}

void sflow_cweno_flux::jqmin(lexer *p,fdm2D *b, slice& f, int ipol)
{   
    q1 = f(i,j-2);
	q2 = f(i,j-1);
	q3 = f(i,j);
	q4 = f(i,j+1);
	q5 = f(i,j+2);
}

void sflow_cweno_flux::iqmax(lexer *p,fdm2D *b, slice& f, int ipol)
{
    q1 = f(i+3,j);
	q2 = f(i+2,j);
	q3 = f(i+1,j);
	q4 = f(i,j);
	q5 = f(i-1,j);
}

void sflow_cweno_flux::jqmax(lexer *p,fdm2D *b, slice& f, int ipol)
{
	q1 = f(i,j+3);
	q2 = f(i,j+2);
	q3 = f(i,j+1);
	q4 = f(i,j);
	q5 = f(i,j-1);
}

void sflow_cweno_flux::is_min(fdm2D *b, slice& f)
{
	is1 = tttw*pow(q1 - 2.0*q2 + q3, 2.0) + fourth*pow(q1 - 4.0*q2 + 3.0*q3, 2.0);
	is2 = tttw*pow(q2 - 2.0*q3 + q4, 2.0) + fourth*pow(q2 - q4, 2.0);
	is3 = tttw*pow(q3 - 2.0*q4 + q5, 2.0) + fourth*pow(3.0*q3 - 4.0*q4 + q5, 2.0);
}

void sflow_cweno_flux::is_max(fdm2D *b, slice& f)
{
	is1 = tttw*pow(q1 - 2.0*q2 + q3, 2.0) + fourth*pow(q1 - 4.0*q2 + 3.0*q3, 2.0);
	is2 = tttw*pow(q2 - 2.0*q3 + q4, 2.0) + fourth*pow(q2 - q4, 2.0);
	is3 = tttw*pow(q3 - 2.0*q4 + q5, 2.0) + fourth*pow(3.0*q3 - 4.0*q4 + q5, 2.0);
}

void sflow_cweno_flux::alpha()
{
	alpha1=tenth/pow(epsilon+is1,2.0);
	alpha2=sixten/pow(epsilon+is2,2.0);
	alpha3=treten/pow(epsilon+is3,2.0);
}

void sflow_cweno_flux::weight()
{
	w1=alpha1/(alpha1+alpha2+alpha3);
	w2=alpha2/(alpha1+alpha2+alpha3);
	w3=alpha3/(alpha1+alpha2+alpha3);
}
