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

#include"geo_raycast.h"
#include"lexer.h"
#include<algorithm>

// Vertical rays through the column centres (FNPF), taken from the FNPF 6DOF body
// geometry (sixdof_obj::ray_cast_fnpf): the crossings with the trimesh are an exact
// inside/outside test (parity) for every node of the column.

void geo_raycast::column_hits(lexer *p, double **tri_x, double **tri_y, double **tri_z, int ts, int te,
                              double dxm, vector<vector<double> > &hits)
{
    int i,j,n;
    
    const int ny = p->knoy;
    
    hits.assign(p->knox*p->knoy,vector<double>());
    
    // 2D: the surface is extruded in y, cast in its mid-plane
    double ymid=0.0;
    
    if(p->j_dir==0 && te>ts)
    {
        double ylo=1.0e20, yhi=-1.0e20;
        
        for(n=ts; n<te; ++n)
        for(int q=0; q<3; ++q)
        {
        ylo = MIN(ylo,tri_y[n][q]);
        yhi = MAX(yhi,tri_y[n][q]);
        }
        
        ymid = 0.5*(ylo+yhi);
    }
    
    // small irrational offset against rays through edges and vertices
    const double ex = 0.1234567891e-6*dxm;
    const double ey = 0.3141592653e-6*dxm;
    
    for(n=ts; n<te; ++n)
    {
        const double x0 = tri_x[n][0], x1 = tri_x[n][1], x2 = tri_x[n][2];
        const double y0 = tri_y[n][0], y1 = tri_y[n][1], y2 = tri_y[n][2];
        const double z0 = tri_z[n][0], z1 = tri_z[n][1], z2 = tri_z[n][2];
        
        const double det = (y1-y2)*(x0-x2) + (x2-x1)*(y0-y2);
        
        // vertical facets are never crossed by a vertical ray
        if(fabs(det)<1.0e-30)
        continue;
        
        const double txmin = MIN(x0,MIN(x1,x2));
        const double txmax = MAX(x0,MAX(x1,x2));
        const double tymin = MIN(y0,MIN(y1,y2));
        const double tymax = MAX(y0,MAX(y1,y2));
        
        ILOOP
        {
            const double xr = p->XP[IP] + ex;
            
            if(xr<txmin || xr>txmax)
            continue;
            
            JLOOP
            {
                const double yr = ((p->j_dir==0) ? ymid : p->YP[JP]) + ey;
                
                if(yr<tymin || yr>tymax)
                continue;
                
                const double l0 = ((y1-y2)*(xr-x2) + (x2-x1)*(yr-y2))/det;
                const double l1 = ((y2-y0)*(xr-x2) + (x0-x2)*(yr-y2))/det;
                const double l2 = 1.0 - l0 - l1;
                
                if(l0<0.0 || l1<0.0 || l2<0.0)
                continue;
                
                hits[i*ny+j].push_back(l0*z0 + l1*z1 + l2*z2);
            }
        }
    }
    
    // a ray through a shared edge or vertex hits every adjacent facet:
    // count coincident crossings once
    const double tol = 1.0e-9*dxm;
    
    for(vector<double> &h : hits)
    {
        if(h.empty())
        continue;
        
        sort(h.begin(),h.end());
        
        size_t nu=1;
        
        for(size_t m=1; m<h.size(); ++m)
        if(h[m]-h[nu-1]>tol)
        h[nu++] = h[m];
        
        h.resize(nu);
    }
}
