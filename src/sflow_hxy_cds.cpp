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

#include"sflow_hxy_cds.h"
#include"lexer.h"
#include"fdm2D.h"
#include"slice.h"
#include"sflow_flux_face_FOU.h"
#include"sflow_flux_face_CDS.h"
#include"sflow_flux_HJ_CDS.h"
#include"patchBC_interface.h"

sflow_hxy_cds::sflow_hxy_cds(lexer* p, patchBC_interface *ppBC)  
{
    pBC = ppBC;
    
    if(p->A216==1)
    pflux = new sflow_flux_face_FOU(p);
        
    if(p->A216==2)
    pflux = new sflow_flux_face_CDS(p);
    
}

sflow_hxy_cds::~sflow_hxy_cds()
{
}

void sflow_hxy_cds::start(lexer* p, slice& hx, slice& hy, slice& depth, int *wet, slice& eta, slice& uvel, slice& vvel)
{
	double eps=0.0;

    SLICELOOP1
	{
	pflux->u_flux(4,uvel,ivel1,ivel2);
/*
	if(ivel1>eps)
    hx(i,j) = eta(i,j) +  MIN(depth(i,j), depth(i+1,j));
	
	if(ivel1<-eps)
    hx(i,j) = eta(i+1,j) +  MIN(depth(i,j), depth(i+1,j));
	
	if(fabs(ivel1)<=eps)*/
    hx(i,j) = 0.5*(eta(i,j)+eta(i+1,j)) + 0.5*(depth(i,j)+depth(i+1,j));
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
	
	/*if(jvel1>eps)
    hy(i,j) = eta(i,j) + MIN(depth(i,j), depth(i,j+1));
	
	if(jvel1<-eps)
    hy(i,j) = eta(i,j+1) + MIN(depth(i,j), depth(i,j+1));
	
	if(fabs(jvel1)<=eps)*/
    hy(i,j) = 0.5*(eta(i,j)+eta(i,j+1)) + 0.5*(depth(i,j)+depth(i,j+1));
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


