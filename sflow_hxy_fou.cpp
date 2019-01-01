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

#include"sflow_hxy_fou.h"
#include"lexer.h"
#include"slice.h"
#include"sflow_flux_face_FOU.h"
#include"sflow_flux_face_CDS.h"
#include"sflow_flux_HJ_CDS.h"

sflow_hxy_fou::sflow_hxy_fou(lexer* p) 
{
    if(p->A216==1)
    pflux = new sflow_flux_face_FOU(p);
        
    if(p->A216>=2)
    pflux = new sflow_flux_face_CDS(p);
}

sflow_hxy_fou::~sflow_hxy_fou()
{
}

void sflow_hxy_fou::start(lexer* p, slice& hx, slice& hy, slice& depth, slice& eta, slice& uvel, slice& vvel)
{
	double eps=1.0e-7;
	
    SLICELOOP1
	{
	pflux->u_flux(4,uvel,ivel1,ivel2);

	if(ivel1>eps)
    hx(i,j) = eta(i,j) + depth(i,j);
	
	if(ivel1<-eps)
    hx(i,j) = eta(i+1,j) + depth(i+1,j);
	
	if(fabs(ivel1)<=eps)
    hx(i,j) = MAX(eta(i,j),eta(i+1,j)) + MIN(depth(i,j), depth(i+1,j));
	}
	
	
	SLICELOOP2
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


