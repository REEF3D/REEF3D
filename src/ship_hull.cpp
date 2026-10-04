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

#include"ship_hull.h"
#include<cmath>

int ship_hull::clip_below(const double *a, const double *b, const double *c, double zw, double (*out)[3])
{
    // Sutherland-Hodgman against the half space z <= zw
    const double *v[3] = {a,b,c};
    int n=0;
    
    for(int q=0; q<3; ++q)
    {
        const double *s = v[q];
        const double *e = v[(q+1)%3];
        
        const bool sin = s[2]<=zw;
        const bool ein = e[2]<=zw;
        
        if(sin)
        {
            out[n][0]=s[0]; out[n][1]=s[1]; out[n][2]=s[2];
            ++n;
        }
        
        if(sin!=ein)
        {
            const double t = (zw - s[2])/(e[2] - s[2]);
            out[n][0] = s[0] + t*(e[0]-s[0]);
            out[n][1] = s[1] + t*(e[1]-s[1]);
            out[n][2] = zw;
            ++n;
        }
    }
    
    return n;
}

double ship_hull::polygon_area(const double (*v)[3], int n)
{
    if(n<3)
    return 0.0;
    
    // half the norm of the summed cross products of a fan from vertex 0
    double sx=0.0, sy=0.0, sz=0.0;
    
    for(int q=1; q<n-1; ++q)
    {
        const double ux = v[q][0]-v[0][0], uy = v[q][1]-v[0][1], uz = v[q][2]-v[0][2];
        const double wx = v[q+1][0]-v[0][0], wy = v[q+1][1]-v[0][1], wz = v[q+1][2]-v[0][2];
        
        sx += uy*wz - uz*wy;
        sy += uz*wx - ux*wz;
        sz += ux*wy - uy*wx;
    }
    
    return 0.5*sqrt(sx*sx + sy*sy + sz*sz);
}

double ship_hull::wetted_surface(double **tx, double **ty, double **tz, int tricount, double zw, bool skip_y)
{
    double S=0.0;
    double poly[4][3];
    
    for(int n=0; n<tricount; ++n)
    {
        const double a[3] = {tx[n][0],ty[n][0],tz[n][0]};
        const double b[3] = {tx[n][1],ty[n][1],tz[n][1]};
        const double c[3] = {tx[n][2],ty[n][2],tz[n][2]};
        
        if(skip_y)
        {
            const double ux=b[0]-a[0], uy=b[1]-a[1], uz=b[2]-a[2];
            const double wx=c[0]-a[0], wy=c[1]-a[1], wz=c[2]-a[2];
            const double nx=uy*wz-uz*wy, ny=uz*wx-ux*wz, nz=ux*wy-uy*wx;
            const double nn=sqrt(nx*nx+ny*ny+nz*nz);
            
            if(nn>0.0 && fabs(ny)>0.9*nn)
            continue;
        }
        
        const int np = clip_below(a,b,c,zw,poly);
        S += polygon_area(poly,np);
    }
    
    return S;
}

void ship_hull::waterline_extent(double **tx, double **ty, double **tz, int tricount, double zw, double &xa, double &xf)
{
    double poly[4][3];
    
    xa =  1.0e20;
    xf = -1.0e20;
    
    for(int n=0; n<tricount; ++n)
    {
        const double a[3] = {tx[n][0],ty[n][0],tz[n][0]};
        const double b[3] = {tx[n][1],ty[n][1],tz[n][1]};
        const double c[3] = {tx[n][2],ty[n][2],tz[n][2]};
        
        const int np = clip_below(a,b,c,zw,poly);
        
        for(int q=0; q<np; ++q)
        {
            xa = poly[q][0]<xa ? poly[q][0] : xa;
            xf = poly[q][0]>xf ? poly[q][0] : xf;
        }
    }
    
    if(xa>xf)
    xa = xf = 0.0;
}

void ship_hull::draft_strips(double **tx, double **ty, double **tz, int tricount, double zw, double xa, double xf, int nstrip,
                             std::vector<double> &xs, std::vector<double> &dx, std::vector<double> &T)
{
    xs.assign(nstrip,0.0);
    dx.assign(nstrip,0.0);
    T.assign(nstrip,0.0);
    
    if(nstrip<1 || xf<=xa)
    return;
    
    const double h = (xf-xa)/double(nstrip);
    std::vector<double> zmin(nstrip,zw);
    double poly[4][3];
    
    for(int n=0; n<tricount; ++n)
    {
        const double a[3] = {tx[n][0],ty[n][0],tz[n][0]};
        const double b[3] = {tx[n][1],ty[n][1],tz[n][1]};
        const double c[3] = {tx[n][2],ty[n][2],tz[n][2]};
        
        const int np = clip_below(a,b,c,zw,poly);
        
        if(np<3)
        continue;
        
        double pxmin=1.0e20, pxmax=-1.0e20, pzmin=1.0e20;
        
        for(int q=0; q<np; ++q)
        {
            pxmin = poly[q][0]<pxmin ? poly[q][0] : pxmin;
            pxmax = poly[q][0]>pxmax ? poly[q][0] : pxmax;
            pzmin = poly[q][2]<pzmin ? poly[q][2] : pzmin;
        }
        
        // strips overlapped by the polygon
        int is = int(floor((pxmin-xa)/h));
        int ie = int(floor((pxmax-xa)/h));
        is = is<0 ? 0 : (is>nstrip-1 ? nstrip-1 : is);
        ie = ie<0 ? 0 : (ie>nstrip-1 ? nstrip-1 : ie);
        
        for(int s=is; s<=ie; ++s)
        zmin[s] = pzmin<zmin[s] ? pzmin : zmin[s];
    }
    
    for(int s=0; s<nstrip; ++s)
    {
        xs[s] = xa + (double(s)+0.5)*h;
        dx[s] = h;
        T[s]  = zw - zmin[s];
    }
}
