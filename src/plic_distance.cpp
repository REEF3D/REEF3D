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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#include"plic_distance.h"
#include"lexer.h"
#include"ghostcell.h"
#include"increment.h"
#include<cmath>
#include<algorithm>
#include<vector>

namespace
{
    const double unset = 1.0e20;
}

plic_distance::plic_distance(lexer *p) : dist(p)
{
}

double plic_distance::plane_constant(double nx, double ny, double nz, double dx, double dy, double dz, double V0)
{
    // Scardovelli & Zaleski (2000), analytic inverse; the same as VOF_PLIC::reconstructPlane_alt
    double n_a = std::fabs(nx)*dx;
    double n_b = std::fabs(ny)*dy;
    double n_c = std::fabs(nz)*dz;
    const double s2 = std::sqrt(n_a*n_a + n_b*n_b + n_c*n_c);
    if(s2<=0.0)
    return 0.0;
    n_a/=s2; n_b/=s2; n_c/=s2;

    double n_1,n_2,n_3;
    if(n_b>=n_a) { n_1=n_a; n_2=n_b; }
    else         { n_1=n_b; n_2=n_a; }

    if(n_c>=n_2)
    n_3=n_c;
    else if(n_c<n_1)
    { n_3=n_2; n_2=n_1; n_1=n_c; }
    else
    { n_3=n_2; n_2=n_c; }

    const double V = 0.5-std::fabs(V0-0.5);
    double r;

    if(n_1+n_2<=2.0*V*n_3)                                                             // case 5
    r = V*n_3 + 0.5*(n_1+n_2);

    else if(3.0*n_2*(2.0*V*n_3+n_1-n_2)<=n_1*n_1 && n_1*n_1<=6.0*V*n_2*n_3)              // case 2
    r = 0.5*n_1 + std::sqrt(2.0*V*n_2*n_3 - n_1*n_1/12.0);

    else if(6.0*V*n_2*n_3<n_1*n_1)                                                      // case 1
    r = std::cbrt(6.0*V*n_1*n_2*n_3);

    else                                                                                // cases 3, 4
    {
        const double x3  = 81.0*n_1*n_2*(n_1+n_2-2.0*n_3*V);
        const double y32 = std::fdim(23328.0*n_1*n_1*n_1*n_2*n_2*n_2, x3*x3);
        const double u3  = std::cbrt(x3*x3+y32);
        const double f3  = n_1+n_2-(7.5595264*n_1*n_2+0.26456684*u3)/std::sqrt(u3)*std::sin(0.5235988-std::atan(std::sqrt(y32)/x3)/3.0);

        if(f3<=n_3)
        r = f3;
        else
        {
            const double t4  = 9.0*(n_1+n_2+n_3)*(n_1+n_2+n_3)-18.0;
            const double x4  = std::max(n_1*n_2*n_3*(324.0-648.0*V),1.0e-35);
            const double y42 = std::fdim(4.0*t4*t4*t4, x4*x4);
            const double u4  = std::cbrt(x4*x4+y42);
            r = 0.5*(n_1+n_2+n_3)-(0.20998684*t4+0.13228342*u4)/std::sqrt(u4)*std::sin(0.5235988-std::atan(std::sqrt(y42)/x4)/3.0);
        }
    }

    const double r0 = (V0-0.5>=0.0) ? 0.5*(n_1+n_2+n_3)-r : -(0.5*(n_1+n_2+n_3)-r);

    return r0*s2;
}

namespace
{
    struct pvec { double x,y,z; };
    inline pvec operator-(pvec a, pvec b) { return {a.x-b.x, a.y-b.y, a.z-b.z}; }
    inline pvec operator+(pvec a, pvec b) { return {a.x+b.x, a.y+b.y, a.z+b.z}; }
    inline pvec operator*(double s, pvec a) { return {s*a.x, s*a.y, s*a.z}; }
    inline double dot(pvec a, pvec b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
    inline pvec cross(pvec a, pvec b) { return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x}; }

    double segment_distance(pvec p, pvec a, pvec b)
    {
        const pvec ab = b-a;
        const double l2 = dot(ab,ab);
        const double t = l2>0.0 ? std::min(std::max(dot(p-a,ab)/l2,0.0),1.0) : 0.0;
        const pvec d = p-(a+t*ab);
        return std::sqrt(dot(d,d));
    }

    // the polygon where the plane n.(x - c) = alpha cuts the cell (centre c, sizes h), ordered
    // around its centroid; empty when the plane misses the cell
    std::vector<pvec> facet(pvec n, double alpha, pvec c, pvec h)
    {
        std::vector<pvec> v;
        const double sx[2]={-0.5,0.5};
        auto f = [&](pvec q){return dot(n,q-c)-alpha;};

        for(int axis=0; axis<3; ++axis)
        for(int a=0; a<2; ++a)
        for(int b=0; b<2; ++b)
        {
            pvec q0, q1;
            if(axis==0) { q0={c.x-0.5*h.x, c.y+sx[a]*h.y, c.z+sx[b]*h.z}; q1={c.x+0.5*h.x, q0.y, q0.z}; }
            if(axis==1) { q0={c.x+sx[a]*h.x, c.y-0.5*h.y, c.z+sx[b]*h.z}; q1={q0.x, c.y+0.5*h.y, q0.z}; }
            if(axis==2) { q0={c.x+sx[a]*h.x, c.y+sx[b]*h.y, c.z-0.5*h.z}; q1={q0.x, q0.y, c.z+0.5*h.z}; }
            const double f0=f(q0), f1=f(q1);
            if((f0<=0.0 && f1>0.0) || (f0>0.0 && f1<=0.0))
            v.push_back(q0+(f0/(f0-f1))*(q1-q0));
        }
        if(v.size()<2)
        return {};

        pvec g={0.0,0.0,0.0};
        for(const pvec &q : v)
        g = g+q;
        g = (1.0/v.size())*g;
        const pvec t = std::fabs(n.x)<0.9 ? pvec{1.0,0.0,0.0} : pvec{0.0,1.0,0.0};
        pvec u = cross(n,t);
        u = (1.0/std::sqrt(dot(u,u)))*u;
        const pvec w = cross(n,u);
        std::sort(v.begin(),v.end(),[&](const pvec &a, const pvec &b)
        {return std::atan2(dot(a-g,w),dot(a-g,u)) < std::atan2(dot(b-g,w),dot(b-g,u));});
        return v;
    }

    // the distance of p to the convex polygon v in the plane of normal n
    double facet_distance(pvec p, pvec n, const std::vector<pvec> &v)
    {
        if(v.size()==2)
        return segment_distance(p,v[0],v[1]);

        const double height = dot(n,p-v[0]);
        const pvec q = p-height*n;
        int pos=0, neg=0;
        for(size_t e=0; e<v.size(); ++e)
        {
            const double s = dot(cross(v[(e+1)%v.size()]-v[e],q-v[e]),n);
            if(s>1.0e-14) ++pos;
            if(s<-1.0e-14) ++neg;
        }
        if(pos==0 || neg==0)
        return std::fabs(height);

        double d=1.0e20;
        for(size_t e=0; e<v.size(); ++e)
        d = std::min(d,segment_distance(p,v[e],v[(e+1)%v.size()]));
        return d;
    }
}

void plic_distance::candidate(lexer *p, int ii, int jj, int kk, double nx, double ny, double nz, double px, double py, double pz)
{
    // the facet of the plane n.(x - P) = 0 (P on it) in cell (ii,jj,kk): its distance to the
    // cells of the band around it. A plane through a face of the cell (P on the face, n its
    // normal): the face itself. In 2D the cells of the row (y) are one slab: the facet is the
    // segment in the x-z plane through the cell centre.
    const int m = increment::marge;
    const int jb = p->j_dir==1 ? band : 0;
    const bool three = p->j_dir==1;

    const pvec n = {nx, three ? ny : 0.0, nz};
    const pvec c = {p->XP[ii+m], three ? p->YP[jj+m] : 0.0, p->ZP[kk+m]};
    const pvec h = {p->DXN[ii+m], three ? p->DYN[jj+m] : 1.0, p->DZN[kk+m]};
    const pvec P = {px, three ? py : 0.0, pz};

    std::vector<pvec> v;
    const double ax=std::fabs(n.x), ay=std::fabs(n.y), az=std::fabs(n.z);
    const bool face = (ax==1.0 || ay==1.0 || az==1.0) && std::fabs(std::fabs(dot(n,P-c))-0.5*(ax*h.x+ay*h.y+az*h.z))<1.0e-12*(h.x+h.y+h.z);

    if(face)
    {
        // the four corners of the face, in order
        const pvec a = ax==1.0 ? pvec{0.0,0.5*h.y,0.0} : pvec{0.5*h.x,0.0,0.0};
        const pvec b = az==1.0 ? pvec{0.0,0.5*h.y,0.0} : pvec{0.0,0.0,0.5*h.z};
        v = {P-a-b, P+a-b, P+a+b, P-a+b};
    }
    else
    v = facet(n, dot(n,P-c), c, h);

    if(!three)
    {
        // the segment: the points of the polygon at y = 0 (two of its corners each)
        std::vector<pvec> s;
        for(const pvec &q : v)
        {
            const pvec r = {q.x, 0.0, q.z};
            bool known=false;
            for(const pvec &o : s)
            if(std::fabs(o.x-r.x)+std::fabs(o.z-r.z)<1.0e-12*(h.x+h.z))
            known=true;
            if(!known)
            s.push_back(r);
        }
        v = s;
    }
    if(v.size()<2)
    return;

    for(int q=ii-band; q<=ii+band; ++q)
    for(int r=jj-jb; r<=jj+jb; ++r)
    for(int s=kk-band; s<=kk+band; ++s)
    if(q>=0 && q<p->knox && r>=0 && r<p->knoy && s>=0 && s<p->knoz)
    {
        const pvec x = {p->XP[q+m], three ? p->YP[r+m] : 0.0, p->ZP[s+m]};
        const double d = facet_distance(x,n,v);

        if(d<dist(q,r,s))
        dist(q,r,s) = d;
    }
}

field4& plic_distance::start(lexer *p, ghostcell *pgc, field &f)
{
    const int m = increment::marge;
    const bool three = p->j_dir==1;
    const int jlo = three ? -2 : 0;
    const int jhi = three ? p->knoy+1 : p->knoy-1;

    auto fluid = [p](int q, int r, int s)
    {
        return p->flag4[(q-p->imin)*p->jmax*p->kmax + (r-p->jmin)*p->kmax + s-p->kmin]>0;
    };

    for(int q=0; q<p->knox; ++q)
    for(int r=0; r<p->knoy; ++r)
    for(int s=0; s<p->knoz; ++s)
    dist(q,r,s) = unset;

    const double w[3] = {1.0,2.0,1.0};

    // the planes of the interface cells, the ghost cells of the rank too (their planes reach
    // the cells of the band in the rank)
    for(int q=-2; q<=p->knox+1; ++q)
    for(int r=jlo; r<=jhi; ++r)
    for(int s=-2; s<=p->knoz+1; ++s)
    {
        if(!fluid(q,r,s))
        continue;

        const double v = f(q,r,s);
        const double xc = p->XP[q+m], yc = p->YP[r+m], zc = p->ZP[s+m];

        if(v>=p->F93 && v<=p->F94)
        {
            // Youngs: weighted central differences over the 3x3(x3) stencil, n = -grad(vof)
            double gx=0.0, gy=0.0, gz=0.0, sw=0.0;
            const int jn = three ? 1 : 0;

            for(int b=-jn; b<=jn; ++b)
            for(int c=-1; c<=1; ++c)
            {
                const double wt = (three ? w[b+1] : 1.0)*w[c+1];
                gx += wt*(f(q+1,r+b,s+c)-f(q-1,r+b,s+c));
                gz += wt*(f(q+c,r+b,s+1)-f(q+c,r+b,s-1));
                if(three)
                gy += wt*(f(q+c,r+1,s+b)-f(q+c,r-1,s+b));
                sw += wt;
            }
            gx /= sw*(p->XP[q+1+m]-p->XP[q-1+m]);
            gz /= sw*(p->ZP[s+1+m]-p->ZP[s-1+m]);
            if(three)
            gy /= sw*(p->YP[r+1+m]-p->YP[r-1+m]);

            const double g = std::sqrt(gx*gx+gy*gy+gz*gz);
            if(g<1.0e-12)
            continue;

            const double nx=-gx/g, ny=-gy/g, nz=-gz/g;
            const double a = plane_constant(nx,ny,nz,p->DXN[q+m],three?p->DYN[r+m]:1.0,p->DZN[s+m],v);

            // a point of the plane: the centre moved by alpha along n
            candidate(p,q,r,s,nx,ny,nz,xc+a*nx,yc+a*ny,zc+a*nz);
        }

        // a full cell next to an empty one: the plane of their face
        if(v>p->F94)
        {
            const int dq[6] = {-1,1,0,0,0,0}, dr[6] = {0,0,-1,1,0,0}, ds[6] = {0,0,0,0,-1,1};

            for(int e=0; e<6; ++e)
            {
                if(!three && dr[e]!=0)
                continue;

                const int q2=q+dq[e], r2=r+dr[e], s2=s+ds[e];
                if(!fluid(q2,r2,s2) || f(q2,r2,s2)>=p->F93)
                continue;

                const double px = dq[e]==0 ? xc : 0.5*(xc+p->XP[q2+m]);
                const double py = dr[e]==0 ? yc : 0.5*(yc+p->YP[r2+m]);
                const double pz = ds[e]==0 ? zc : 0.5*(zc+p->ZP[s2+m]);
                candidate(p,q,r,s,double(dq[e]),double(dr[e]),double(ds[e]),px,py,pz);
            }
        }
    }

    // the sign of each cell's own vof; outside the band the vertical distance to the free
    // surface of its column, at least the width of the band
    for(int q=0; q<p->knox; ++q)
    for(int r=0; r<p->knoy; ++r)
    {
        double ztop=-1.0e20, air=0.0;
        for(int s=0; s<p->knoz; ++s)
        if(fluid(q,r,s))
        {
            ztop = std::max(ztop, p->ZP[s+m]+0.5*p->DZN[s+m]);
            air += (1.0-f(q,r,s))*p->DZN[s+m];
        }
        const double zs = ztop-air;

        for(int s=0; s<p->knoz; ++s)
        {
            const double sgn = f(q,r,s)>=0.5 ? 1.0 : -1.0;

            if(dist(q,r,s)>=unset)
            {
                const double far = (band+0.5)*std::max(p->DXN[q+m],p->DZN[s+m]);
                dist(q,r,s) = sgn*std::max(std::fabs(zs-p->ZP[s+m]),far);
            }
            else
            dist(q,r,s) = sgn*dist(q,r,s);
        }
    }

    pgc->start4(p,dist,1);

    return dist;
}
