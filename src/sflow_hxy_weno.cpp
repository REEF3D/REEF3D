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

#include"sflow_hxy_weno.h"
#include"lexer.h"
#include"fdm2D.h"
#include"slice.h"
#include"sflow_flux_face_FOU.h"
#include"sflow_flux_face_CDS.h"
#include"sflow_flux_face_HJ.h"
#include"patchBC_interface.h"

sflow_hxy_weno::sflow_hxy_weno(lexer* p, patchBC_interface *ppBC) :tttw(13.0/12.0),fourth(1.0/4.0),third(1.0/3.0),
			sevsix(7.0/6.0),elvsix(11.0/6.0),sixth(1.0/6.0),fivsix(5.0/6.0),tenth(1.0/10.0),
			sixten(6.0/10.0),treten(3.0/10.0),epsilon(0.000001),smallnum(1.0e-20)
{
    pBC = ppBC;
    
    if(p->A216==1)
    pflux = new sflow_flux_face_FOU(p);
        
    if(p->A216==2)
    pflux = new sflow_flux_face_CDS(p);
    
    if(p->A216==4)
    pflux = new sflow_flux_face_HJ(p);
}

sflow_hxy_weno::~sflow_hxy_weno()
{
}

void sflow_hxy_weno::start(lexer* p, slice& hx, slice& hy, slice& depth, int *wet, slice& eta, slice& uvel, slice& vvel)
{
    double eps=1.0e-7;

	SLICELOOP1
	{
	pflux->u_flux(4,uvel,ivel1,ivel2);
    if(ivel1>=0.0)
    hx(i,j) = fx(p,eta,1,ivel1) + MIN(depth(i,j), depth(i+1,j));
    
    if(ivel1<0.0)
    hx(i,j) = fx(p,eta,1,ivel1) + MIN(depth(i,j), depth(i+1,j));
	}
    
    
    if(p->F50==1 || p->F50==4)
    for(n=0;n<p->gcslout_count;n++)
    {
    i=p->gcslout[n][0];
    j=p->gcslout[n][1];
    
        if(wet[IJ]==1)
        {
        pflux->u_flux(4,uvel,ivel1,ivel2);

        if(ivel1>eps)
        hx(i,j) = eta(i,j) + depth(i,j);
        
        if(ivel1<-eps)
        hx(i,j) = eta(i+1,j) + depth(i+1,j);
        
        if(fabs(ivel1)<=eps)
        hx(i,j) = MAX(eta(i,j),eta(i+1,j)) + MIN(depth(i,j), depth(i+1,j));
        }
    }
    
    int qq;    
    for(qq=0;qq<pBC->obj_count;++qq)
    if(pBC->patch[qq]->waterlevel_flag==0)
    for(n=0;n<pBC->patch[qq]->gcb_count;++n)
    if(pBC->patch[qq]->gcb[n][3]==1 || pBC->patch[qq]->gcb[n][3]==4)
    {
    if(pBC->patch[qq]->gcb[n][3]==1)
    i=pBC->patch[qq]->gcb[n][0]-1;
    
    j=pBC->patch[qq]->gcb[n][1];

        
        if(wet[IJ]==1)
        {
        pflux->u_flux(4,uvel,ivel1,ivel2);

        if(ivel1>eps)
        hx(i,j) = eta(i,j) + depth(i,j);
        
        if(ivel1<-eps)
        hx(i,j) = eta(i+1,j) + depth(i+1,j);
        
        if(fabs(ivel1)<=eps)
        hx(i,j) = MAX(eta(i,j),eta(i+1,j)) + MIN(depth(i,j), depth(i+1,j));
        }
    }
	
	SLICELOOP2
	{
	pflux->v_flux(4,vvel,jvel1,jvel2);
    if(jvel1>=0.0)
	hy(i,j) = fy(p,eta,2,jvel1) + MIN(depth(i,j), depth(i,j+1));
    
    if(jvel1<0.0)
	hy(i,j) = fy(p,eta,2,jvel1) + MIN(depth(i,j), depth(i,j+1));
    }

    for(qq=0;qq<pBC->obj_count;++qq)
    if(pBC->patch[qq]->waterlevel_flag==0)
    for(n=0;n<pBC->patch[qq]->gcb_count;++n)
    if(pBC->patch[qq]->gcb[n][3]==3 || pBC->patch[qq]->gcb[n][3]==2)
    {
    
    i=pBC->patch[qq]->gcb[n][0];
    
    if(pBC->patch[qq]->gcb[n][3]==3)
    j=pBC->patch[qq]->gcb[n][1]-1;

        
        if(wet[IJ]==1)
        {
        pflux->v_flux(4,vvel,jvel1,jvel2);
	
        if(jvel1>eps)
        hy(i,j) = eta(i,j) + depth(i,j);
        
        if(jvel1<-eps)
        hy(i,j) = eta(i,j+1) + depth(i,j+1);
        
        if(fabs(jvel1)<=eps)
        hy(i,j) = MAX(eta(i,j),eta(i,j+1)) + MIN(depth(i,j), depth(i,j+1));
        }
    }
	 
}

double sflow_hxy_weno::fx(lexer *p, slice& f, int ipol, double advec)
{
    grad = 0.0;

	if(advec>=0.0)
	{
	iqmin(p,f,ipol);
	is(f);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	if(advec<0.0)
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

double sflow_hxy_weno::fy(lexer *p, slice& f, int ipol, double advec)
{
    grad = 0.0;

	if(advec>=0.0)
	{
	jqmin(p,f,ipol);
	is(f);
	alpha();
	weight();
	
	grad = (w1*( q1*third - q2*sevsix + q3*elvsix)
	      + w2*(-q2*sixth + q3*fivsix + q4*third)
		  + w3*( q3*third + q4*fivsix - q5*sixth));
	}

	if(advec<0.0)
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


void sflow_hxy_weno::iqmin(lexer *p, slice& f, int ipol)
{	
	q1 = f(i-2,j);
	q2 = f(i-1,j);
	q3 = f(i,j);
	q4 = f(i+1,j);
	q5 = f(i+2,j);
}

void sflow_hxy_weno::jqmin(lexer *p, slice& f, int ipol)
{
	q1 = f(i,j-2);
	q2 = f(i,j-1);
	q3 = f(i,j);
	q4 = f(i,j+1);
	q5 = f(i,j+2);
}

void sflow_hxy_weno::iqmax(lexer *p, slice& f, int ipol)
{
	q1 = f(i+3,j);
	q2 = f(i+2,j);
	q3 = f(i+1,j);
	q4 = f(i,j);
	q5 = f(i-1,j);
}

void sflow_hxy_weno::jqmax(lexer *p, slice& f, int ipol)
{
	q1 = f(i,j+3);
	q2 = f(i,j+2);
	q3 = f(i,j+1);
	q4 = f(i,j);
	q5 = f(i,j-1);
}

void sflow_hxy_weno::is(slice& f)
{
	is1 = tttw*pow(q1 - 2.0*q2 + q3, 2.0) + fourth*pow(q1 - 4.0*q2 + 3.0*q3, 2.0);
	is2 = tttw*pow(q2 - 2.0*q3 + q4, 2.0) + fourth*pow(q2 - q4, 2.0);
	is3 = tttw*pow(q3 - 2.0*q4 + q5, 2.0) + fourth*pow(3.0*q3 - 4.0*q4 + q5, 2.0);
}

void sflow_hxy_weno::alpha()
{
	alpha1=tenth/pow(epsilon+is1,2.0);
	alpha2=sixten/pow(epsilon+is2,2.0);
	alpha3=treten/pow(epsilon+is3,2.0);
}

void sflow_hxy_weno::weight()
{
	w1=alpha1/(alpha1+alpha2+alpha3);
	w2=alpha2/(alpha1+alpha2+alpha3);
	w3=alpha3/(alpha1+alpha2+alpha3);
}
