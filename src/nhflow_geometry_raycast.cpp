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

#include"nhflow_geometry.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#define WLVL (fabs(d->WL(i,j))>0.00005?d->WL(i,j):1.0e20)

void nhflow_geometry::ray_cast(lexer *p, fdm_nhf *d, ghostcell *pgc, double *LS)
{
    // a patch of the mesh refinement: the bounds of the ray cast are the patch arrays
    // (amr_ray_box), the extent of the lexer restored at the end
    struct bounds
    {
        lexer *p;
        bool on;
        double ox,ex,oy,ey;
        bounds(lexer *pp, bool o, const double *b) : p(pp), on(o)
        {
            if(!on)
            return;
            ox=p->originx; ex=p->endx; oy=p->originy; ey=p->endy;
            p->originx=b[0]; p->endx=b[1]; p->originy=b[2]; p->endy=b[3];
        }
        ~bounds()
        {
            if(on)
            {
            p->originx=ox; p->endx=ex; p->originy=oy; p->endy=ey;
            }
        }
    } bnd(p,amr_grid,amr_bx);
    
    zmin = 1.0e8;
    zmax = -1.0e8;
    

    LOOP
    WETDRY
    {
    zmin = MIN(zmin, p->ZSP[IJK]);
    zmax = MAX(zmax, p->ZSP[IJK]);
    }
    
    LOOP
	{
    IO[IJK]=1;
	LS[IJK]=1.0e8;
	}
    
    	
    // geometry core kernels (geo_raycast)
    for(int rayiter=0; rayiter<2; ++rayiter)
    {
        for(int qn=0;qn<entity_sum;++qn)
        {
            if(rayiter==0)
            {
            const int mode = (qn<int(ent_raymode.size())) ? ent_raymode[qn] : 1;
            
            georay.sigma_io(p,tri_x,tri_y,tri_z,tstart[qn],tend[qn],DSM,mode,IO,CL,CR);
            
            if(qn<int(ent_invert.size()) && ent_invert[qn]==1)
            LOOP
            IO[IJK] = -IO[IJK];
            }
            
            if(rayiter==1)
            {
            pgc->startintV(p,IO,1);
            
            // VRANS: distance to the exposed surface only (layer_faces), all faces otherwise
            georay.sigma_band(p,tri_x,tri_y,tri_z,tstart[qn],tend[qn],NB*DSM,LS,
                              tri_layer.empty() ? nullptr : tri_layer.data());
            }
        }
    }
    
    // cell centres on the surface (distance at round-off level) get LS = 0 exactly: the parity
    // test perturbs the ray in a fixed direction, so such cells would be inside on one side of
    // the solid and outside on the other (asymmetric level set of a symmetric solid)
    const double ls_eps = 1.0e-10*DSM;
    
    LOOP
    WETDRY
    {
        if(IO[IJK]==-1)
        LS[IJK]=-fabs(LS[IJK]);
        
        
        if(IO[IJK]==1)
        LS[IJK]=fabs(LS[IJK]);
        
        if(fabs(LS[IJK])<ls_eps)
        LS[IJK]=0.0;
    }
	
	LOOP
    WETDRY
	{
		if(LS[IJK]>100.0*DSM)
		LS[IJK]=100.0*DSM;
		
		if(LS[IJK]<-100.0*DSM)
		LS[IJK]=-100.0*DSM;
	}
    
    LOOP
    if(p->wet[IJ]==0)
    LS[IJK]=100.0*DSM;
    
    
	pgc->start5V(p,LS,1); 
}

void nhflow_geometry::entity_box(int n, double *b) const
{
    b[0]=b[2]=1.0e20;
    b[1]=b[3]=-1.0e20;
    
    for(int q=tstart[n]; q<tend[n]; ++q)
    for(int v=0; v<3; ++v)
    {
        b[0] = MIN(b[0],tri_x[q][v]);
        b[1] = MAX(b[1],tri_x[q][v]);
        b[2] = MIN(b[2],tri_y[q][v]);
        b[3] = MAX(b[3],tri_y[q][v]);
    }
}
