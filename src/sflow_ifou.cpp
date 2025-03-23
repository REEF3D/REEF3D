/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"sflow_ifou.h"
#include"lexer.h"
#include"fdm2D.h"
#include"sflow_flux_face_FOU.h"
#include"sflow_flux_face_CDS.h"
#include"sflow_flux_face_HJ.h"
#include"sflow_flux_HJ_CDS.h"

sflow_ifou::sflow_ifou(lexer* p)
{
    if(p->A216==1)
    pflux = new sflow_flux_face_FOU(p);
        
    if(p->A216==2)
    pflux = new sflow_flux_face_CDS(p);
    
    if(p->A216==4)
    pflux = new sflow_flux_face_HJ(p);
}

sflow_ifou::~sflow_ifou()
{
}

void sflow_ifou::start(lexer* p, fdm2D* b, slice& f, int ipol, slice& uvel, slice& vvel)
{
    count=0;
    
    if(ipol==1)
    SLICELOOP1
    aij(p,b,f,1,uvel,vvel);

    if(ipol==2)
    SLICELOOP2
    aij(p,b,f,2,uvel,vvel);
    
    if(ipol==4)
    SLICELOOP4
    aij(p,b,f,4,uvel,vvel);
}

void sflow_ifou::aij(lexer* p,fdm2D* b,slice& f,int ipol, slice& uvel, slice& vvel)
{
    udir=vdir=0.0;
    
    pflux->u_flux(ipol,uvel,iadvec,ivel2);
    pflux->v_flux(ipol,vvel,jadvec,jvel2);


	if(0.5*(ivel1+ivel2)>=0.0)
    udir=1.0;
    
    if(0.5*(jvel1+jvel2)>=0.0)
    vdir=1.0;
    

	 
	 b->M.p[count] = udir*ivel2/p->DXM - (1.0-udir)*ivel1/p->DXM
					+ (vdir*jvel2/p->DXM - (1.0-vdir)*jvel1/p->DXM)*p->y_dir;

	 
	 b->M.s[count] = -udir*ivel1/p->DXM;
	 b->M.n[count] =  (1.0-udir)*ivel2/p->DXM;
	 
	 b->M.e[count] = -vdir*jvel1/p->DXM*p->y_dir;
	 b->M.w[count] =  (1.0-vdir)*jvel2/p->DXM*p->y_dir;
	 
	 ++count;
}

