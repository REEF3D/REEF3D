/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"6DOF_obj_cfd.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"fieldint.h"

void sixdof_obj_cfd::ray_cast(lexer *p, fdm *a, ghostcell *pgc)
{
	LOOP
	{
    fbio(i,j,k)=1;
	a->fb(i,j,k)=1.0e8;
	}
	
    // geometry core kernels (geo_raycast), subdomain interior
    const geo_cart g = geo_raycast::interior(p);
    
    for(rayiter=0; rayiter<2; ++rayiter)
    {

        for(int qn=0;qn<entity_sum;++qn)
        {
            if(rayiter==0)
            {
            georay.cart_io(p,0,tri_x,tri_y,tri_z,tstart[qn],tend[qn],g,fbio.V,cutl.V,cutr.V);
            
            if(p->j_dir==1)
            georay.cart_io(p,1,tri_x,tri_y,tri_z,tstart[qn],tend[qn],g,fbio.V,cutl.V,cutr.V);
            
            georay.cart_io(p,2,tri_x,tri_y,tri_z,tstart[qn],tend[qn],g,fbio.V,cutl.V,cutr.V);
            }
        
            if(rayiter==1 && p->X188==1)
            {
            pgc->gcparaxint(p,fbio,4);
            
            georay.cart_dist(p,0,tri_x,tri_y,tri_z,tstart[qn],tend[qn],g,fbio.V,a->fb.V);
            
            if(p->j_dir==1)
            georay.cart_dist(p,1,tri_x,tri_y,tri_z,tstart[qn],tend[qn],g,fbio.V,a->fb.V);
            
            georay.cart_dist(p,2,tri_x,tri_y,tri_z,tstart[qn],tend[qn],g,fbio.V,a->fb.V);
            }
            
            if(rayiter==1 && p->X188==2)
            {
            pgc->gcparaxint(p,fbio,4);
            
            georay.cart_vertexdist(p,tri_x,tri_y,tri_z,tstart[qn],tend[qn],a->fb.V);
            }
        }
    }
    
    // cell centres on the surface (distance at round-off level) get fb = 0 exactly: the parity
    // test perturbs the ray in a fixed direction, so such cells would be inside on one side of
    // the body and outside on the other (asymmetric level set of a symmetric body)
    const double fb_eps = 1.0e-10*p->DXM;
    
    LOOP
    {
        if(fbio(i,j,k)==-1)
        a->fb(i,j,k)=-fabs(a->fb(i,j,k));
        
        
        if(fbio(i,j,k)==1)
        a->fb(i,j,k)=fabs(a->fb(i,j,k));
        
        if(fabs(a->fb(i,j,k))<fb_eps)
        a->fb(i,j,k)=0.0;
    }
    
	
	LOOP
	{
		if(a->fb(i,j,k)>10.0*p->DXM)
		a->fb(i,j,k)=10.0*p->DXM;
		
		if(a->fb(i,j,k)<-10.0*p->DXM)
		a->fb(i,j,k)=-10.0*p->DXM;
	}
    
	pgc->start4a(p,a->fb,50); 
}





