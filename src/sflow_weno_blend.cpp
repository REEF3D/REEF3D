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

#include"sflow_weno_blend.h"
#include"lexer.h"
#include"fdm2D.h"
#include"sflow_flux_face_FOU.h"
#include"sflow_flux_face_CDS.h"
#include"sflow_flux_face_HJ.h"

sflow_weno_blend::sflow_weno_blend(lexer* p):tttw(13.0/12.0),fourth(1.0/4.0),third(1.0/3.0),
			sevsix(7.0/6.0),elvsix(11.0/6.0),sixth(1.0/6.0),fivsix(5.0/6.0),tenth(1.0/10.0),
			sixten(6.0/10.0),treten(3.0/10.0),epsilon(0.000001),smallnum(1.0e-20)
{
    if(p->A216==1)
    pflux = new sflow_flux_face_FOU(p);
        
    if(p->A216==2)
    pflux = new sflow_flux_face_CDS(p);
    
    if(p->A216==4)
    pflux = new sflow_flux_face_HJ(p);
    
    pflux_hj = new sflow_flux_face_HJ(p);

}

sflow_weno_blend::~sflow_weno_blend()
{
}

void sflow_weno_blend::start(lexer* p, fdm2D* b, slice& f, int ipol, slice& uvel, slice& vvel)
{
    if(ipol==1)
    SLICELOOP1
    {
    if(p->flagslice1[Im1J]>0 && p->flagslice1[Ip1J]>0)
    b->F(i,j)+=aij_flux(p,b,f,1,uvel,vvel);
    
    if(p->flagslice1[Im1J]<0 || p->flagslice1[Ip1J]<0)
    b->F(i,j)+=aij_hj(p,b,f,1,uvel,vvel);
    }

    if(ipol==2)
    SLICELOOP2
    {
    if(p->flagslice2[IJm1]>0 && p->flagslice2[IJp1]>0)
    b->G(i,j)+=aij_flux(p,b,f,2,uvel,vvel);
    
    if(p->flagslice2[IJm1]<0 || p->flagslice2[IJp1]<0)
    b->G(i,j)+=aij_hj(p,b,f,2,uvel,vvel);
    }
    
    if(ipol==4)
    SLICELOOP4
    b->L(i,j)+=aij_flux(p,b,f,4,uvel,vvel);

    if(ipol==5)
    SLICELOOP4
    b->L(i,j)+=aij_flux(p,b,f,5,uvel,vvel);

}

double sflow_weno_blend::aij_flux(lexer* p,fdm2D* b,slice& f,int ipol, slice& uvel, slice& vvel)
{
		pflux->u_flux(ipol,uvel,ivel1,ivel2);
        pflux->v_flux(ipol,vvel,jvel1,jvel2);

		i-=1;
		fu1 = fx_flux(p,b,f,ipol,ivel1);
		i+=1;
		
		fu2 = fx_flux(p,b,f,ipol,ivel2);

		
		j-=1;
		fv1 = fy_flux(p,b,f,ipol,jvel1);
		j+=1;
		
		fv2 = fy_flux(p,b,f,ipol,jvel2);
		
		
		L =   - ((ivel2*fu2-ivel1*fu1)/p->DXM) 
		      - ((jvel2*fv2-jvel1*fv1)/p->DXM);
  			  
		return L;
}

double sflow_weno_blend::aij_hj(lexer* p,fdm2D* b,slice& f,int ipol, slice& uvel, slice& vvel)
{
		pflux_hj->u_flux(ipol,uvel,iadvec,ivel2);
		pflux_hj->v_flux(ipol,vvel,jadvec,jvel2);

		L = -iadvec*fx_hj(p,b,f,ipol,iadvec) - jadvec*fy_hj(p,b,f,ipol,jadvec);

		return L;
}

double sflow_weno_blend::fx_flux(lexer *p,fdm2D *b, slice& f, int ipol, double advec)
{
    grad = 0.0;

	if(advec>0.0)
	{
	iqmin_flux(p,b,f,ipol);
	is(f);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	if(advec<0.0)
	{
	iqmax_flux(p,b,f,ipol);
	is(f);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	return grad;
}

double sflow_weno_blend::fx_hj(lexer *p,fdm2D *b, slice& f, int ipol, double advec)
{
    grad = 0.0;

	if(advec>0.0)
	{
	iqmin_hj(p,b,f,ipol);
	is(f);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	if(advec<0.0)
	{
	iqmax_hj(p,b,f,ipol);
	is(f);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	return grad;
}

double sflow_weno_blend::fy_flux(lexer *p,fdm2D *b, slice& f, int ipol, double advec)
{
    grad = 0.0;

	if(advec>0.0)
	{
	jqmin_flux(p,b,f,ipol);
	is(f);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	if(advec<0.0)
	{
	jqmax_flux(p,b,f,ipol);
	is(f);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}
	
	return grad;
}

double sflow_weno_blend::fy_hj(lexer *p,fdm2D *b, slice& f, int ipol, double advec)
{
    grad = 0.0;

	if(advec>0.0)
	{
	jqmin_hj(p,b,f,ipol);
	is(f);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	if(advec<0.0)
	{
	jqmax_hj(p,b,f,ipol);
	is(f);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}
	
	return grad;
}


void sflow_weno_blend::iqmin_flux(lexer *p,fdm2D *b, slice& f, int ipol)
{	
	q1 = f(i-2,j);
	q2 = f(i-1,j);
	q3 = f(i,j);
	q4 = f(i+1,j);
	q5 = f(i+2,j);
}

void sflow_weno_blend::jqmin_flux(lexer *p,fdm2D *b, slice& f, int ipol)
{
	q1 = f(i,j-2);
	q2 = f(i,j-1);
	q3 = f(i,j);
	q4 = f(i,j+1);
	q5 = f(i,j+2);
}

void sflow_weno_blend::iqmax_flux(lexer *p,fdm2D *b, slice& f, int ipol)
{
	q1 = f(i+3,j);
	q2 = f(i+2,j);
	q3 = f(i+1,j);
	q4 = f(i,j);
	q5 = f(i-1,j);
}

void sflow_weno_blend::jqmax_flux(lexer *p,fdm2D *b, slice& f, int ipol)
{
	q1 = f(i,j+3);
	q2 = f(i,j+2);
	q3 = f(i,j+1);
	q4 = f(i,j);
	q5 = f(i,j-1);
}

void sflow_weno_blend::iqmin_hj(lexer *p,fdm2D *b, slice& f, int ipol)
{	
	q1 = (f(i-2,j) - f(i-3,j))/p->DXM;
	q2 = (f(i-1,j) - f(i-2,j))/p->DXM;
	q3 = (f(i,j)   - f(i-1,j))/p->DXM;
	q4 = (f(i+1,j) - f(i,j)  )/p->DXM;
	q5 = (f(i+2,j) - f(i+1,j))/p->DXM;
}

void sflow_weno_blend::jqmin_hj(lexer *p,fdm2D *b, slice& f, int ipol)
{
	q1 = (f(i,j-2) - f(i,j-3))/p->DXM;
	q2 = (f(i,j-1) - f(i,j-2))/p->DXM;
	q3 = (f(i,j)   - f(i,j-1))/p->DXM;
	q4 = (f(i,j+1) - f(i,j)  )/p->DXM;
	q5 = (f(i,j+2) - f(i,j+1))/p->DXM;
}

void sflow_weno_blend::iqmax_hj(lexer *p,fdm2D *b, slice& f, int ipol)
{
	q1 = (f(i+3,j) - f(i+2,j))/p->DXM;
	q2 = (f(i+2,j) - f(i+1,j))/p->DXM;
	q3 = (f(i+1,j) - f(i,j)  )/p->DXM;
	q4 = (f(i,j)   - f(i-1,j))/p->DXM;
	q5 = (f(i-1,j) - f(i-2,j))/p->DXM;
}

void sflow_weno_blend::jqmax_hj(lexer *p,fdm2D *b, slice& f, int ipol)
{
	q1 = (f(i,j+3) - f(i,j+2))/p->DXM;
	q2 = (f(i,j+2) - f(i,j+1))/p->DXM;
	q3 = (f(i,j+1) - f(i,j)  )/p->DXM;
	q4 = (f(i,j)   - f(i,j-1))/p->DXM;
	q5 = (f(i,j-1) - f(i,j-2))/p->DXM;
}

void sflow_weno_blend::is(slice& f)
{
	is1 = tttw*pow(q1 - 2.0*q2 + q3, 2.0) + fourth*pow(q1 - 4.0*q2 + 3.0*q3, 2.0);
	is2 = tttw*pow(q2 - 2.0*q3 + q4, 2.0) + fourth*pow(q2 - q4, 2.0);
	is3 = tttw*pow(q3 - 2.0*q4 + q5, 2.0) + fourth*pow(3.0*q3 - 4.0*q4 + q5, 2.0);
}

void sflow_weno_blend::alpha()
{
	alpha1=tenth/pow(epsilon+is1,2.0);
	alpha2=sixten/pow(epsilon+is2,2.0);
	alpha3=treten/pow(epsilon+is3,2.0);
}

void sflow_weno_blend::weight()
{
	w1=alpha1/(alpha1+alpha2+alpha3);
	w2=alpha2/(alpha1+alpha2+alpha3);
	w3=alpha3/(alpha1+alpha2+alpha3);
}
