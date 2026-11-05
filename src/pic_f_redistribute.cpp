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

#include"pic_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include<algorithm>
#include<cmath>
#include<vector>

// The redistribution of the particles (F 78: every F 78 steps; F 81 1, 2), after Kelly et al. 2015,
// section 3.5 (Ando et al. 2012): the particles drift into rows and sheets with gaps between, the
// FLIP particles in particular. Each particle is pushed away from its neighbours within
// re = 2 s (s the spacing of the seeding, h / F 83),
//     dx_i = - F79 s sum_j k(|x_j - x_i|, re) (x_j - x_i) / |x_j - x_i|,   k = (1 - r^2/re^2)^3,
// at most s/4 per iteration, in 4 iterations. All particles move at once from the positions of the
// iteration before (those of the neighbouring subdomains as the halo), so the result does not
// depend on the decomposition; the walls of the domain repel with the mirror images of the
// particles. Fixed are the particles within one cell of the surface (phi < h, as in Kelly et al.,
// so the surface stays where the particles put it) and moves into a solid or a floating body are
// not done. With F 81 1, 2 a particle keeps its velocity plus the change of the water velocity
// between the two positions (no PIC smoothing of the velocity).

void pic_f::redistribute(lexer *p, fdm *a, ghostcell *pgc)
{
    const int iters = 4;
    const double beta = p->F79;
    const double kre = 2.0;
    const bool d3 = p->j_dir==1;
    const int fm = fmarg;
    const int ns = std::max(1,int(std::ceil(kre/double(npd))));

    const int nx = p->knox + 2*fm;
    const int ny = d3 ? p->knoy + 2*fm : 1;
    const int nz = p->knoz + 2*fm;

    auto bin = [&](double x, double y, double z)
    {
        const int ci = cell(p->XN,p->knox,x) + fm;
        const int cj = d3 ? cell(p->YN,p->knoy,y) + fm : 0;
        const int ck = cell(p->ZN,p->knoz,z) + fm;
        return (ci*ny + cj)*nz + ck;
    };

    // the particles that may move: deeper than a cell under the surface
    const size_t np = P.size();
    std::vector<char> mobile(np,0);
    std::vector<double> x0(3*np);

    for(size_t m=0; m<np; ++m)
    {
        const pic_data &q = P[m];
        const int ci = std::min(std::max(cell(p->XN,p->knox,q.x),0),p->knox-1);
        const int cj = d3 ? std::min(std::max(cell(p->YN,p->knoy,q.y),0),p->knoy-1) : 0;
        const int ck = std::min(std::max(cell(p->ZN,p->knoz,q.z),0),p->knoz-1);
        const double h = crad[lid(p,ci,cj,ck)]/2.0;

        mobile[m] = ipol(p,a->phi,4,q.x,q.y,q.z) >= h && q.x>=xglo_s;
        x0[3*m] = q.x;
        x0[3*m+1] = q.y;
        x0[3*m+2] = q.z;
    }

    std::vector<double> Q, H, dx(3*np);
    std::vector<int> head, next;

    for(int it=0; it<iters; ++it)
    {
        halo(p,pgc,H);

        // the neighbours: own, halo and the mirror images across the walls of the domain
        Q.clear();

        auto add = [&](double x, double y, double z)
        {
            Q.push_back(x);
            Q.push_back(y);
            Q.push_back(z);
        };

        // self: the particle itself is in the list already, only its images are added
        auto images = [&](double x, double y, double z, bool self)
        {
            const double bx = 3.0*p->DXM, by = 3.0*p->DXM, bz = 3.0*p->DXM;
            double xs[2] = {x,x}, ys[2] = {y,y}, zs[2] = {z,z};
            int mx=1, my=1, mz=1;

            if(x-xglo_s<bx)
            xs[mx++] = 2.0*xglo_s - x;
            else if(xglo_e-x<bx)
            xs[mx++] = 2.0*xglo_e - x;

            if(d3)
            {
            if(y-yglo_s<by)
            ys[my++] = 2.0*yglo_s - y;
            else if(yglo_e-y<by)
            ys[my++] = 2.0*yglo_e - y;
            }

            if(z-zglo_s<bz)
            zs[mz++] = 2.0*zglo_s - z;
            else if(zglo_e-z<bz)
            zs[mz++] = 2.0*zglo_e - z;

            for(int a1=0; a1<mx; ++a1)
            for(int a2=0; a2<my; ++a2)
            for(int a3=0; a3<mz; ++a3)
            if(!self || a1+a2+a3>0)
            add(xs[a1],ys[a2],zs[a3]);
        };

        // the own particles first: Q[3 m] is particle m
        for(const auto &q : P)
        add(q.x,q.y,q.z);

        for(const auto &q : P)
        images(q.x,q.y,q.z,true);

        for(size_t m=0; m+2<H.size(); m+=3)
        images(H[m],H[m+1],H[m+2],false);

        const int nq = int(Q.size()/3);
        head.assign(size_t(nx)*ny*nz,-1);
        next.assign(nq,-1);

        for(int n=0; n<nq; ++n)
        {
            const int b = bin(Q[3*n],Q[3*n+1],Q[3*n+2]);
            next[n] = head[b];
            head[b] = n;
        }

        for(size_t m=0; m<np; ++m)
        {
            dx[3*m] = dx[3*m+1] = dx[3*m+2] = 0.0;

            if(!mobile[m])
            continue;

            const pic_data &q = P[m];
            const int ci = cell(p->XN,p->knox,q.x);
            const int cj = d3 ? cell(p->YN,p->knoy,q.y) : 0;
            const int ck = cell(p->ZN,p->knoz,q.z);
            const int cic = std::min(std::max(ci,0),p->knox-1);
            const int cjc = std::min(std::max(cj,0),p->knoy-1);
            const int ckc = std::min(std::max(ck,0),p->knoz-1);
            const double s = crad[lid(p,cic,cjc,ckc)]/2.0/double(npd);
            const double re = kre*s;
            const double re2 = re*re;
            double f[3] = {0.0,0.0,0.0};

            for(int bi=std::max(ci-ns,-fm); bi<=std::min(ci+ns,p->knox+fm-1); ++bi)
            for(int bj=(d3 ? std::max(cj-ns,-fm) : 0); bj<=(d3 ? std::min(cj+ns,p->knoy+fm-1) : 0); ++bj)
            for(int bk=std::max(ck-ns,-fm); bk<=std::min(ck+ns,p->knoz+fm-1); ++bk)
            for(int n=head[((bi+fm)*ny + (d3 ? bj+fm : 0))*nz + bk+fm]; n>=0; n=next[n])
            {
                if(n==int(m))
                continue;

                const double rx = Q[3*n]-q.x;
                const double ry = d3 ? Q[3*n+1]-q.y : 0.0;
                const double rz = Q[3*n+2]-q.z;
                const double r2 = rx*rx + ry*ry + rz*rz;

                if(r2>=re2 || r2<1.0e-24*re2)
                continue;

                const double r = std::sqrt(r2);
                const double w = kernel(r2,re2)/r;

                f[0] += w*rx;
                f[1] += w*ry;
                f[2] += w*rz;
            }

            double d[3] = {-beta*s*f[0], -beta*s*f[1], -beta*s*f[2]};
            const double dn = std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
            const double dmax = 0.25*s;

            if(dn>dmax)
            for(int c=0; c<3; ++c)
            d[c] *= dmax/dn;

            dx[3*m] = d[0];
            dx[3*m+1] = d3 ? d[1] : 0.0;
            dx[3*m+2] = d[2];
        }

        // all at once; not into solids, bodies or out of the domain
        for(size_t m=0; m<np; ++m)
        {
            if(dx[3*m]==0.0 && dx[3*m+1]==0.0 && dx[3*m+2]==0.0)
            continue;

            pic_data &q = P[m];
            const double x = q.x + dx[3*m], y = q.y + dx[3*m+1], z = q.z + dx[3*m+2];

            if(x<=xglo_s || x>=xglo_e || z<=zglo_s || z>=zglo_e || (d3 && (y<=yglo_s || y>=yglo_e)))
            continue;

            if(blocked(p,a,x,y,z))
            continue;

            q.x = x;
            q.y = y;
            q.z = z;
        }
    }

    // F 81 1, 2: the change of the water velocity between the old and the new position
    for(size_t m=0; m<np; ++m)
    {
        pic_data &q = P[m];

        if(!flipmode || (q.x==x0[3*m] && q.y==x0[3*m+1] && q.z==x0[3*m+2]))
        continue;

        q.u += ipol(p,*uw,1,q.x,q.y,q.z) - ipol(p,*uw,1,x0[3*m],x0[3*m+1],x0[3*m+2]);
        if(d3)
        q.v += ipol(p,*vw,2,q.x,q.y,q.z) - ipol(p,*vw,2,x0[3*m],x0[3*m+1],x0[3*m+2]);
        q.w += ipol(p,*ww,3,q.x,q.y,q.z) - ipol(p,*ww,3,x0[3*m],x0[3*m+1],x0[3*m+2]);
    }

    xchange(p,pgc);
}
