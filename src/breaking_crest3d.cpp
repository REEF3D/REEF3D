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

#include"breaking_crest3d.h"
#include"lexer.h"
#include"ghostcell.h"
#include"slice.h"
#include<cmath>
#include<algorithm>
#include<sys/stat.h>
#include<mpi.h>
#include"reef3d_mpmd.h"

namespace
{
    // least-squares slope of v over t, the last n of m samples (oldest first)
    double bc3_slope(const double *t, const double *v, int m, int n)
    {
        n = std::min(m,n);

        if(n<2)
        return 0.0;

        const int s = m-n;
        double tm=0.0, vm=0.0;

        for(int q=s; q<m; ++q)
        {
            tm += t[q];
            vm += v[q];
        }
        tm/=double(n);
        vm/=double(n);

        double num=0.0, den=0.0;

        for(int q=s; q<m; ++q)
        {
            num += (t[q]-tm)*(v[q]-vm);
            den += (t[q]-tm)*(t[q]-tm);
        }

        return (den>0.0) ? num/den : 0.0;
    }

    // sin^2 ramp over the first and last 10 % of [0,1]
    double bc3_ramp(double s)
    {
        const double a = 0.1;

        if(s<0.0 || s>1.0)
        return 0.0;

        if(s<a)
        return pow(sin(0.5*M_PI*s/a),2.0);

        if(s>1.0-a)
        return pow(sin(0.5*M_PI*(1.0-s)/a),2.0);

        return 1.0;
    }

    // append to a ring of 8 (oldest first)
    void bc3_push(double *a, double *b, double *c, int &n, double va, double vb, double vc)
    {
        if(n==8)
        {
            for(int q=0; q<7; ++q)
            {
                a[q]=a[q+1];
                b[q]=b[q+1];
                if(c!=nullptr)
                c[q]=c[q+1];
            }
            n=7;
        }

        a[n]=va;
        b[n]=vb;
        if(c!=nullptr)
        c[n]=vc;
        ++n;
    }

    // all points of all ranks
    void bc3_allgather(const std::vector<bc3_point> &loc, std::vector<bc3_point> &all, int nproc)
    {
        const int sz = (int)sizeof(bc3_point);
        int nb = (int)loc.size()*sz;
        std::vector<int> cnt(nproc), dsp(nproc);
        MPI_Allgather(&nb,1,MPI_INT,cnt.data(),1,MPI_INT,reef3d_comm());

        int ntot=0;
        for(int q=0; q<nproc; ++q)
        {
            dsp[q]=ntot;
            ntot+=cnt[q];
        }

        all.resize(ntot/sz);
        MPI_Allgatherv(loc.data(),nb,MPI_BYTE,all.data(),cnt.data(),dsp.data(),MPI_BYTE,reef3d_comm());
    }
}

breaking_crest3d::breaking_crest3d(lexer *p)
{
    MPI_Comm_size(reef3d_comm(),&nproc);
    MPI_Comm_rank(reef3d_comm(),&rank);
}

breaking_crest3d::~breaking_crest3d()
{
}

// cell index containing x (cell centres XP, own cells and ghost cells)
int breaking_crest3d::cell_i(lexer *p, double x)
{
    const int lo = -marge+1, hi = p->knox+marge-2;

    if(x<p->XP[lo+marge] || x>p->XP[hi+marge])
    return -1000000000;

    int a=lo, b=hi;

    while(b-a>1)
    {
        const int m=(a+b)/2;

        if(p->XP[m+marge]<=x)
        a=m;
        else
        b=m;
    }

    return a;
}

int breaking_crest3d::cell_j(lexer *p, double y)
{
    if(p->j_dir==0)
    return 0;

    const int lo = -marge+1, hi = p->knoy+marge-2;

    if(y<p->YP[lo+marge] || y>p->YP[hi+marge])
    return -1000000000;

    int a=lo, b=hi;

    while(b-a>1)
    {
        const int m=(a+b)/2;

        if(p->YP[m+marge]<=y)
        a=m;
        else
        b=m;
    }

    return a;
}

// inside the own subdomain (cell faces)
bool breaking_crest3d::own(lexer *p, double x, double y)
{
    if(x<p->XN[marge] || x>=p->XN[p->knox+marge])
    return false;

    if(p->j_dir==1 && (y<p->YN[marge] || y>=p->YN[p->knoy+marge]))
    return false;

    return true;
}

// bilinear (2D grid: linear) value; nan outside the reach of the ghost cells
double breaking_crest3d::ip(lexer *p, slice &f, double x, double y)
{
    const int ii = cell_i(p,x);
    const int jj = cell_j(p,y);

    if(ii<-marge || jj<-marge)
    return NAN;

    const double a = (x-p->XP[ii+marge])/(p->XP[ii+1+marge]-p->XP[ii+marge]);

    if(p->j_dir==0)
    return (1.0-a)*f(ii,0) + a*f(ii+1,0);

    const double b = (y-p->YP[jj+marge])/(p->YP[jj+1+marge]-p->YP[jj+marge]);

    return (1.0-a)*(1.0-b)*f(ii,jj) + a*(1.0-b)*f(ii+1,jj) + (1.0-a)*b*f(ii,jj+1) + a*b*f(ii+1,jj+1);
}

// walls: the sides without inflow / outflow cells (gcbsl4 side codes: 1 x-, 4 x+, 3 y-, 2 y+)
void breaking_crest3d::walls(lexer *p, ghostcell *pgc)
{
    int open[5] = {0,0,0,0,0};

    for(int q=0; q<p->gcslin_count; ++q)
    {
        const int cs = p->gcslin[q][3];
        open[cs==1 ? 1 : cs==4 ? 2 : cs==3 ? 3 : 4] = 1;
    }

    for(int q=0; q<p->gcslout_count; ++q)
    {
        const int cs = p->gcslout[q][3];
        open[cs==1 ? 1 : cs==4 ? 2 : cs==3 ? 3 : 4] = 1;
    }

    for(int e=1; e<=4; ++e)
    wallside[e] = (pgc->globalimax(open[e])==0 && (e<=2 || p->j_dir==1)) ? 1 : 0;

    wallside[0] = 0;
}

// a ray sample beyond a wall is moved back onto the wall (the ray slides along it); true if moved
bool breaking_crest3d::clampxy(lexer *p, double &x, double &y)
{
    const double ex = 1.0e-9*(p->global_xmax-p->global_xmin);
    const double ey = 1.0e-9*(p->global_ymax-p->global_ymin);
    bool moved = false;

    if(wallside[1]==1 && x<p->global_xmin+ex) {x = p->global_xmin+ex; moved = true;}
    if(wallside[2]==1 && x>p->global_xmax-ex) {x = p->global_xmax-ex; moved = true;}
    if(wallside[3]==1 && y<p->global_ymin+ey) {y = p->global_ymin+ey; moved = true;}
    if(wallside[4]==1 && y>p->global_ymax-ey) {y = p->global_ymax-ey; moved = true;}

    return moved;
}

bool breaking_crest3d::wetcell(lexer *p, int ii, int jj)
{
    return p->wet[(ii-p->imin)*p->jmax + (jj-p->jmin)]==1;
}

// distance across one cell along (nx, ny)
double breaking_crest3d::celld(lexer *p, int ii, int jj, double nx, double ny)
{
    const double dx = p->DXN[ii+marge];

    if(p->j_dir==0)
    return dx;

    const double dy = p->DYN[jj+marge];

    return 1.0/sqrt(nx*nx/(dx*dx) + ny*ny/(dy*dy));
}

void breaking_crest3d::detect(lexer *p, ghostcell *pgc, slice &eta, slice &us, slice &vs, slice &WL, int *zone, double t, slice *phi)
{
    const double g = fabs(p->W22);

    prev.swap(pts);
    pts.clear();

    const double dtp = (tprev>-1.0e19) ? t-tprev : 0.0;
    tprev = t;

    // ------------------------------------------------------------
    // 1. crest points
    // ------------------------------------------------------------
    for(i=0; i<p->knox; ++i)
    for(j=0; j<p->knoy; ++j)
    {
        if(p->flagslice4[IJ]<=0 || p->wet[IJ]==0 || (zone!=nullptr && zone[IJ]==1))
        continue;


        int allwet=1;

        for(int qi=-1; qi<=1; ++qi)
        for(int qj=-p->j_dir; qj<=p->j_dir; ++qj)
        if(!wetcell(p,i+qi,j+qj))
        allwet=0;

        if(allwet==0)
        continue;

        // ridge normal
        double nx=1.0, ny=0.0;

        const double x0=p->XP[IP-1], x1=p->XP[IP], x2=p->XP[IP+1];
        const double hxx = 2.0*((eta(i+1,j)-eta(i,j))/(x2-x1) - (eta(i,j)-eta(i-1,j))/(x1-x0))/(x2-x0);

        if(p->j_dir==1)
        {
            const double y0=p->YP[JP-1], y1=p->YP[JP], y2=p->YP[JP+1];
            const double hyy = 2.0*((eta(i,j+1)-eta(i,j))/(y2-y1) - (eta(i,j)-eta(i,j-1))/(y1-y0))/(y2-y0);
            const double hxy = (eta(i+1,j+1)-eta(i+1,j-1)-eta(i-1,j+1)+eta(i-1,j-1))/((x2-x0)*(y2-y0));

            // most negative eigenvalue of [[hxx,hxy],[hxy,hyy]] and its eigenvector
            const double tr = 0.5*(hxx+hyy);
            const double dt = sqrt(0.25*(hxx-hyy)*(hxx-hyy) + hxy*hxy);
            const double lmin = tr-dt;
            const double lmax = tr+dt;

            // a crest: the most negative curvature dominates (ridge, or a peak with lmin ~ lmax < 0); this
            // rejects saddles and the flat ridges across a trough or a front face (noise in the other direction)
            if(lmin>=0.0 || lmin>-0.5*fabs(lmax))
            continue;

            // eigenvector: the larger of (hxy, lmin-hxx) and (lmin-hyy, hxy) (both vanish for a diagonal Hessian)
            const double ax = hxy, ay = lmin-hxx;
            const double bx = lmin-hyy, by = hxy;

            if(ax*ax+ay*ay >= bx*bx+by*by)
            {
                nx = ax;
                ny = ay;
            }
            else
            {
                nx = bx;
                ny = by;
            }

            if(nx*nx+ny*ny<1.0e-300)
            {
                nx = (hxx<=hyy) ? 1.0 : 0.0;
                ny = (hxx<=hyy) ? 0.0 : 1.0;
            }

            const double nn = sqrt(nx*nx+ny*ny);
            nx/=nn;
            ny/=nn;
        }
        else
        if(hxx>=0.0)
        continue;

        // maximum along n: parabola through -d, 0, +d
        const double d = celld(p,i,j,nx,ny);
        const double xc0 = p->XP[IP];
        const double yc0 = (p->j_dir==1) ? p->YP[JP] : 0.0;
        const double e0 = eta(i,j);
        const double em = ip(p,eta,xc0-d*nx,yc0-d*ny);
        const double ep = ip(p,eta,xc0+d*nx,yc0+d*ny);

        if(std::isnan(em) || std::isnan(ep))
        continue;

        // the points along n within the domain (not beyond the outer boundaries, e.g. into the ghost cells of the
        // wave generation)
        {
            const double gx0 = p->global_xmin, gx1 = p->global_xmax;
            const double xm = xc0-d*nx, xp = xc0+d*nx;

            if(xm<gx0 || xp<gx0 || xm>gx1 || xp>gx1)
            continue;

            if(p->j_dir==1)
            {
                const double gy0 = p->global_ymin, gy1 = p->global_ymax;
                const double ym = yc0-d*ny, yp = yc0+d*ny;

                if(ym<gy0-1.0e-12 || yp<gy0-1.0e-12 || ym>gy1+1.0e-12 || yp>gy1+1.0e-12)
                continue;
            }
        }

        if(!((e0>em && e0>=ep) || (e0>=em && e0>ep)))
        continue;

        const double cv = em-2.0*e0+ep;
        double sc = 0.0, ec = e0;

        if(cv<-1.0e-14)
        {
            sc = 0.5*d*(em-ep)/cv;
            sc = MAX(sc,-0.5*d);
            sc = MIN(sc, 0.5*d);
            ec = e0 - (em-ep)*(em-ep)/(8.0*cv);
        }

        bc3_point P;
        P.x = xc0 + sc*nx;
        P.y = yc0 + sc*ny;
        P.nx = nx;
        P.ny = ny;
        P.eta = ec;
        P.wl = MAX(WL(i,j),0.0);
        P.d = d;
        P.i = i;
        P.j = j;
        P.rank = rank;

        // continue the nearest crest point of the last step within two cells (same orientation)
        const double rm = 2.0*MAX(p->DXN[IP], (p->j_dir==1) ? p->DYN[JP] : 0.0);
        const bc3_point *best=nullptr;
        double bd = rm*rm;

        for(int list=0; list<2; ++list)
        for(const bc3_point &Q : (list==0 ? prev : halo))
        {
            // predicted position of Q
            const double mv = (Q.hasC==1) ? Q.C*dtp : 0.0;
            const double dx = Q.x+mv*Q.nx-P.x;
            const double dy = Q.y+mv*Q.ny-P.y;
            const double dd = dx*dx+dy*dy;

            if(dd<=bd && fabs(Q.nx*nx+Q.ny*ny)>0.7)
            {
                bd = dd;
                best = &Q;
            }
        }

        if(best!=nullptr)
        {
            const bc3_point &Q = *best;

            // orientation of the crest point it continues
            if(Q.nx*P.nx+Q.ny*P.ny<0.0)
            {
                P.nx=-P.nx;
                P.ny=-P.ny;
            }

            P.id = Q.id;
            P.state = Q.state;
            P.ton = Q.ton;
            P.b = Q.b;
            P.C = Q.C;
            P.hasC = Q.hasC;
            P.nh = Q.nh;
            P.nbh = Q.nbh;

            for(int q=0; q<8; ++q)
            {
                P.th[q]=Q.th[q];
                P.xh[q]=Q.xh[q];
                P.yh[q]=Q.yh[q];
                P.tb[q]=Q.tb[q];
                P.bh[q]=Q.bh[q];
            }
        }
        else
        {
            P.id = (long long)rank*1000000000LL + nextid;
            ++nextid;
        }

        // on the front face of a breaking crest point of the last step: continue its breaking
        if(P.state==0)
        {
            const bc3_point *bq=nullptr;
            double bs = 1.0e20;
            const double wl = 1.5*MAX(p->DXN[IP], (p->j_dir==1) ? p->DYN[JP] : 0.0);

            for(int list=0; list<2; ++list)
            for(const bc3_point &Q : (list==0 ? prev : halo))
            if(Q.state==1 && Q.sa>0.0 && fabs(Q.nx*nx+Q.ny*ny)>0.7)
            {
                const double mv = (Q.hasC==1) ? Q.C*dtp : 0.0;
                const double dx = P.x-(Q.x+mv*Q.nx);
                const double dy = P.y-(Q.y+mv*Q.ny);
                const double s = dx*Q.nx + dy*Q.ny;
                const double q = fabs(-dx*Q.ny + dy*Q.nx);

                if(s>=-Q.d && s<=Q.sa+Q.d && q<=wl && fabs(s)<bs)
                {
                    bs = fabs(s);
                    bq = &Q;
                }
            }

            if(bq!=nullptr)
            {
                if(bq->nx*P.nx+bq->ny*P.ny<0.0)
                {
                    P.nx=-P.nx;
                    P.ny=-P.ny;
                }

                P.state = 1;
                P.ton = bq->ton;
                P.b = bq->b;

                if(P.hasC==0)
                {
                    P.C = bq->C;
                    P.hasC = 1;
                }

                best = bq;
            }
        }

        // surface velocity along n at the crest (n as found, before the orientation): the parabola through -d, 0, +d
        // of phi (derivative) or of (us, vs).n, at the crest offset sc
        double usn = 0.0;

        if(phi!=nullptr)
        {
            const double fm = ip(p,*phi,xc0-d*nx,yc0-d*ny);
            const double f0 = (*phi)(i,j);
            const double fp = ip(p,*phi,xc0+d*nx,yc0+d*ny);

            usn = (fp-fm)/(2.0*d) + (fp-2.0*f0+fm)/(d*d)*sc;
        }
        else
        {
            auto un = [&](double x, double y) { return ip(p,us,x,y)*nx + ((p->j_dir==1) ? ip(p,vs,x,y)*ny : 0.0); };
            const double um = un(xc0-d*nx,yc0-d*ny);
            const double u0 = us(i,j)*nx + ((p->j_dir==1) ? vs(i,j)*ny : 0.0);
            const double up = un(xc0+d*nx,yc0+d*ny);

            usn = u0 + (up-um)/(2.0*d)*sc + 0.5*(up-2.0*u0+um)/(d*d)*sc*sc;
        }

        // the ridge normal within 45 deg of the surface velocity: at a travelling crest u_s points along the direction
        // of travel; ridges across the direction of travel (on front faces and in troughs, where the curvature along
        // the direction of travel changes sign) are not crests
        if(p->j_dir==1)
        {
            const double ux = ip(p,us,P.x,P.y), uy = ip(p,vs,P.x,P.y);
            const double um = sqrt(ux*ux+uy*uy);

            if(um>1.0e-12 && fabs(ux*nx+uy*ny)<0.7071*um)
            continue;
        }

        // the sign of usn refers to the n of the detection; P.n may have been turned by the matching
        if(P.nx*nx+P.ny*ny<0.0)
        usn = -usn;

        if(best==nullptr && usn<0.0)
        {
            P.nx=-P.nx;
            P.ny=-P.ny;
            usn=-usn;
        }

        // crest speed: least-squares slope of the positions projected on n; a sample every quarter cell of travel
        const double cs = (P.hasC==1 && fabs(P.C)>0.05) ? fabs(P.C) : MAX(sqrt(g*MAX(P.wl,1.0e-3)), 0.05);
        const double dts = 0.25*d/cs;
        const bool sample = (P.nh==0 || t-P.th[P.nh-1]>=dts);

        if(sample)
        bc3_push(P.th,P.xh,P.yh,P.nh,t,P.x,P.y);

        double tt[9], ss[9];
        int m=0;

        for(int q=0; q<P.nh; ++q)
        {
            tt[m] = P.th[q];
            ss[m] = P.xh[q]*P.nx + P.yh[q]*P.ny;
            ++m;
        }

        if(!sample)
        {
            tt[m] = t;
            ss[m] = P.x*P.nx + P.y*P.ny;
            ++m;
        }

        // C from the own samples; a crest point that took over a breaking crest keeps that C until it has four
        if(m>=4)
        {
            P.hasC = 1;
            P.C = bc3_slope(tt,ss,m,9);
        }
        else
        if(P.state!=1)
        {
            P.hasC = 0;
            P.C = 0.0;
        }

        if(P.hasC==1 && fabs(P.C)<1.0e-6)
        P.hasC=0;

        // a crest travelling against n (strong opposing current): turn n
        if(P.hasC==1 && P.C<0.0)
        {
            P.nx=-P.nx;
            P.ny=-P.ny;
            P.C=-P.C;
            usn=-usn;
        }

        P.us = usn;
        P.B = (P.hasC==1) ? P.us/P.C : 0.0;

        if(P.hasC==1 && sample)
        bc3_push(P.tb,P.bh,nullptr,P.nbh,t,P.B,0.0);

        pts.push_back(P);
    }

    exchange_halo(p);
}

// the crest points within reach of a subdomain border, for the matching of the next step on the other ranks
void breaking_crest3d::exchange_halo(lexer *p)
{
    std::vector<bc3_point> loc;

    for(const bc3_point &P : pts)
    {
        const double r = 3.0*MAX(p->DXN[P.i+marge], (p->j_dir==1) ? p->DYN[P.j+marge] : 0.0);

        if(P.x-r<p->XN[marge] || P.x+r>p->XN[p->knox+marge]
        || (p->j_dir==1 && (P.y-r<p->YN[marge] || P.y+r>p->YN[p->knoy+marge])))
        loc.push_back(P);
    }

    std::vector<bc3_point> all;
    bc3_allgather(loc,all,nproc);

    // keep the ones of the other ranks within reach of the own subdomain
    halo.clear();

    const double xa = p->XN[marge], xb = p->XN[p->knox+marge];
    const double ya = (p->j_dir==1) ? p->YN[marge] : -1.0e20, yb = (p->j_dir==1) ? p->YN[p->knoy+marge] : 1.0e20;

    for(const bc3_point &Q : all)
    if(Q.rank!=rank)
    {
        const double r = 3.0*Q.d+1.0e-12;

        if(Q.x>xa-r && Q.x<xb+r && Q.y>ya-r && Q.y<yb+r)
        halo.push_back(Q);
    }
}

void breaking_crest3d::rays(lexer *p, ghostcell *pgc, slice &eta, slice &WL, double t,
                            const std::function<double(const bc3_point&)> &sback,
                            const std::function<double(const bc3_point&, double, double)> &integrand,
                            const std::function<double(const bc3_point&)> &epsf,
                            slice &vb, int *brk, slice *kw)
{
    const double g = fabs(p->W22);

    if(wallside[0]==-1)
    walls(p,pgc);

    // an axis with a wall side: the ray is not cut there (clampxy), the own subdomain is tested per sample
    const bool xwall = wallside[1]==1 || wallside[2]==1;
    const bool ywall = wallside[3]==1 || wallside[4]==1;

    // ------------------------------------------------------------
    // candidates of all ranks
    // ------------------------------------------------------------
    std::vector<bc3_point> loc;

    for(const bc3_point &P : pts)
    if(P.hasC==1 && (P.B>Bcand || P.state==1))
    loc.push_back(P);

    bc3_allgather(loc,cand,nproc);

    const int nc = (int)cand.size();

    // longest ray: the diagonal of the domain
    const double Lx = p->global_xmax-p->global_xmin;
    const double Ly = (p->j_dir==1) ? p->global_ymax-p->global_ymin : 0.0;
    const double smax_all = sqrt(Lx*Lx+Ly*Ly);

    // parameter interval of a ray in the own subdomain (+ one sample)
    const double xa0 = p->XN[marge], xb0 = p->XN[p->knox+marge];
    const double ya0 = (p->j_dir==1) ? p->YN[marge] : -1.0e20, yb0 = (p->j_dir==1) ? p->YN[p->knoy+marge] : 1.0e20;

    auto range = [&](const bc3_point &R, double dir, double ds, int &m0, int &m1, double mg=0.0)
    {
        // s in [0, smax_all] with x_c + dir s n inside [xa,xb) x [ya,yb), widened by mg
        const double xa = xa0-mg, xb = xb0+mg, ya = ya0-mg, yb = yb0+mg;
        double s0=0.0, s1=smax_all;
        const double vx = dir*R.nx, vy = dir*R.ny;

        if(!xwall)
        {
        if(fabs(vx)>1.0e-14)
        {
            const double a = (xa-R.x)/vx, b = (xb-R.x)/vx;
            s0 = MAX(s0,MIN(a,b));
            s1 = MIN(s1,MAX(a,b));
        }
        else
        if(R.x<xa || R.x>=xb)
        s1=-1.0;
        }

        if(p->j_dir==1 && !ywall)
        {
            if(fabs(vy)>1.0e-14)
            {
                const double a = (ya-R.y)/vy, b = (yb-R.y)/vy;
                s0 = MAX(s0,MIN(a,b));
                s1 = MIN(s1,MAX(a,b));
            }
            else
            if(R.y<ya || R.y>=yb)
            s1=-1.0;
        }

        m0 = (int)ceil(s0/ds - 1.0e-9);
        m1 = (int)floor(s1/ds + 1.0e-9);
    };

    // ------------------------------------------------------------
    // pass 1: first minimum ahead (dir +1) and behind (dir -1); a dry sample ends the ray at the sample before
    // ------------------------------------------------------------
    std::vector<int> ma(MAX(nc,1),1000000000), mb(MAX(nc,1),1000000000);

    for(int k=0; k<nc; ++k)
    {
        const bc3_point &R = cand[k];
        const double ds = 0.5*R.d;

        for(int side=0; side<2; ++side)
        {
            const double dir = (side==0) ? 1.0 : -1.0;
            int m0, m1;
            range(R,dir,ds,m0,m1);
            m0 = MAX(m0,1);

            int found = 1000000000;

            for(int m=m0; m<=m1 && found==1000000000; ++m)
            {
                double x = R.x + dir*m*ds*R.nx;
                double y = R.y + dir*m*ds*R.ny;
                double xm = x-dir*ds*R.nx, ym = y-dir*ds*R.ny;
                double xp = x+dir*ds*R.nx, yp = y+dir*ds*R.ny;
                clampxy(p,x,y);
                clampxy(p,xm,ym);
                clampxy(p,xp,yp);

                if(!own(p,x,y))
                continue;

                const int ii = cell_i(p,x), jj = cell_j(p,y);
                const bool dry = (ii<-marge || jj<-marge || !wetcell(p,ii,jj) || !wetcell(p,ii+1,jj)
                               || (p->j_dir==1 && (!wetcell(p,ii,jj+1) || !wetcell(p,ii+1,jj+1))));

                if(dry)
                {
                    found = m-1;
                    break;
                }

                const double e0 = ip(p,eta,x,y);
                const double em = ip(p,eta,xm,ym);
                const double ep = ip(p,eta,xp,yp);

                if(!std::isnan(ep) && e0<=em && e0<ep)
                found = m;
            }

            if(side==0)
            ma[k]=found;
            else
            mb[k]=found;
        }
    }

    if(nc>0)
    {
        MPI_Allreduce(MPI_IN_PLACE,ma.data(),nc,MPI_INT,MPI_MIN,reef3d_comm());
        MPI_Allreduce(MPI_IN_PLACE,mb.data(),nc,MPI_INT,MPI_MIN,reef3d_comm());
    }

    // painted interval [-s_back, s_ahead]
    for(int k=0; k<nc; ++k)
    {
        bc3_point &R = cand[k];
        const double ds = 0.5*R.d;

        R.sa = (ma[k]<1000000000) ? ma[k]*ds : -1.0;
        R.sb = (mb[k]<1000000000) ? mb[k]*ds : -1.0;
        R.L = (R.sa>0.0 && R.sb>0.0) ? R.sa+R.sb : 2.0*MAX(R.sa,0.0);
        R.kw = (R.L>0.0) ? 2.0*M_PI/R.L : 0.0;
        R.sback = (R.sa>0.0) ? MAX(sback(R),0.0) : 0.0;
    }

    // ------------------------------------------------------------
    // pass 2: trough values, largest front-face slope, dissipation per unit viscosity
    // ------------------------------------------------------------
    std::vector<double> et(MAX(nc,1),0.0), ht(MAX(nc,1),0.0), sl(MAX(nc,1),0.0), I(MAX(nc,1),0.0);

    for(int k=0; k<nc; ++k)
    {
        const bc3_point &R = cand[k];

        if(R.sa<=0.0)
        continue;

        const double ds = 0.5*R.d;
        const int mtr = ma[k];
        const int mback = (int)floor(R.sback/ds + 1.0e-9);

        // ahead (m >= 0) and behind (m < 0)
        for(int side=0; side<2; ++side)
        {
            const double dir = (side==0) ? 1.0 : -1.0;
            int m0, m1;
            range(R,dir,ds,m0,m1);

            if(side==0)
            {
                m0 = MAX(m0,0);
                m1 = MIN(m1,mtr);
            }
            else
            {
                m0 = MAX(m0,1);
                m1 = MIN(m1,mback);
            }

            for(int m=m0; m<=m1; ++m)
            {
                const double s = dir*m*ds;
                double x = R.x + s*R.nx;
                double y = R.y + s*R.ny;
                double xm = x-R.d*R.nx, ym = y-R.d*R.ny;
                double xp = x+R.d*R.nx, yp = y+R.d*R.ny;
                clampxy(p,x,y);
                clampxy(p,xm,ym);
                clampxy(p,xp,yp);

                if(!own(p,x,y))
                continue;

                if(side==0 && m==mtr)
                {
                    et[k] = ip(p,eta,x,y);
                    ht[k] = ip(p,WL,x,y);
                }

                if(side==0)
                {
                    const double em = ip(p,eta,xm,ym);
                    const double ep = ip(p,eta,xp,yp);

                    if(!std::isnan(em) && !std::isnan(ep))
                    sl[k] = MAX(sl[k], fabs(ep-em)/(2.0*R.d));
                }

                const double S = bc3_ramp((s+R.sback)/MAX(R.sa+R.sback,1.0e-12));
                const double w = (m==0 || (side==0 && m==mtr) || (side==1 && m==mback)) ? 0.5 : 1.0;

                if(S>0.0)
                I[k] += w*S*integrand(R,x,y)*ds;
            }
        }
    }

    if(nc>0)
    {
        MPI_Allreduce(MPI_IN_PLACE,et.data(),nc,MPI_DOUBLE,MPI_SUM,reef3d_comm());
        MPI_Allreduce(MPI_IN_PLACE,ht.data(),nc,MPI_DOUBLE,MPI_SUM,reef3d_comm());
        MPI_Allreduce(MPI_IN_PLACE,sl.data(),nc,MPI_DOUBLE,MPI_MAX,reef3d_comm());
        MPI_Allreduce(MPI_IN_PLACE,I.data(),nc,MPI_DOUBLE,MPI_SUM,reef3d_comm());
    }

    // ------------------------------------------------------------
    // state and strength (the same on all ranks)
    // ------------------------------------------------------------
    nbreak = 0;

    for(int k=0; k<nc; ++k)
    {
        bc3_point &R = cand[k];

        R.H = R.eta - et[k];
        R.ht = MAX(ht[k],1.0e-3);
        R.dpn = I[k];

        if(R.sa<=0.0 || R.H<=Hmin)
        {
            R.state=0;
            R.nu=0.0;
            continue;
        }

        // a crest of a wave travels close to the phase speed: maxima that stall (ripples, a flat crest whose maximum
        // jumps) have a small C and a large B. No onset below Cmin of the linear phase speed sqrt(g tanh(kh)/k),
        // k = 2 pi / L, h the depth at the crest.
        // With the period of the generated waves, c_lin is that of the waves at the crest depth: short waves riding on
        // the waves (NHFLOW: L ~ h) have their own crests with u_s of the long wave and a small C.
        double clin = sqrt(g*MAX(R.wl,1.0e-3));

        if(Tdom>0.0)
        {
            const double w = 2.0*M_PI/Tdom;
            const double hh = MAX(R.wl,1.0e-3);
            double kk = MAX(w*w/g, w/sqrt(g*hh));

            for(int q=0; q<50; ++q)
            {
                const double th = tanh(kk*hh);
                const double f = g*kk*th - w*w;
                const double df = g*th + g*kk*hh*(1.0-th*th);
                const double dk = f/df;
                kk -= dk;

                if(fabs(dk)<1.0e-12*kk)
                break;
            }

            clin = w/kk;
        }
        else
        if(R.L>0.0)
        {
            const double kk = 2.0*M_PI/R.L;
            clin = sqrt(g*tanh(kk*MAX(R.wl,1.0e-3))/kk);
        }

        // onset
        if(R.state==0 && R.B>Bon && R.C>=Cmin*clin)
        {
            R.state = 1;
            R.ton = t;
            deepest_onset = MAX(deepest_onset, R.wl-R.eta);

            // deep-water strength from dB/dt at onset (Derakhti, Banner & Kirby 2018)
            R.b = 0.0;

            if(deep==1)
            {
                const double Tb = R.L/MAX(fabs(R.C),1.0e-6);
                const double Gm = Tb*bc3_slope(R.tb,R.bh,R.nbh,5);

                if(Gm>0.30)
                R.b = 0.034*pow(Gm-0.30,2.5);
            }
        }

        if(R.state!=1)
        {
            R.nu=0.0;
            continue;
        }

        // end (not in the step of onset)
        if(R.ton<t && sl[k]<slope_end)
        {
            R.state=0;
            R.nu=0.0;
            continue;
        }

        R.eps = epsf(R);
        R.nu = (R.dpn>1.0e-12) ? R.eps/R.dpn : 0.0;
        R.nu = MIN(R.nu, fabs(R.C)*(R.sa+R.sback));
        R.nu = MAX(R.nu, 0.0);

        ++nbreak;
    }

    // the state goes back to the own crest points
    for(bc3_point &P : pts)
    for(const bc3_point &R : cand)
    if(R.rank==rank && R.i==P.i && R.j==P.j)
    {
        P.state = R.state;
        P.ton = R.ton;
        P.b = R.b;
        P.nu = R.nu;
        P.sa = R.sa;
        break;
    }

    // the halo copies of the points of the other ranks carry the new state as well
    for(bc3_point &Q : halo)
    for(const bc3_point &R : cand)
    if(R.rank==Q.rank && R.i==Q.i && R.j==Q.j)
    {
        Q.state = R.state;
        Q.ton = R.ton;
        Q.b = R.b;
        Q.sa = R.sa;
        break;
    }

    // ------------------------------------------------------------
    // painting
    // ------------------------------------------------------------
    for(int k=0; k<nc; ++k)
    {
        const bc3_point &R = cand[k];

        if(R.state!=1 || R.sa<=0.0)
        continue;

        const double ds = 0.5*R.d;
        const double sf0 = -R.sback-2.0*R.d, sf1 = R.sa+2.0*R.d;
        const double mg = 2.0*R.d + 2.0*MAX(p->DXM, 1.0e-12);

        for(int side=0; side<2; ++side)
        {
        const double dir = (side==0) ? 1.0 : -1.0;
        int m0, m1;
        range(R,dir,ds,m0,m1,mg);
        m0 = MAX(m0, side);
        m1 = MIN(m1, (side==0) ? (int)ceil(sf1/ds) : (int)ceil(-sf0/ds));

        for(int m=m0; m<=m1; ++m)
        {
            const double s = dir*m*ds;
            double x = R.x + s*R.nx;
            double y = R.y + s*R.ny;
            const bool slid = clampxy(p,x,y);

            const int ic = cell_i(p,x), jc = cell_j(p,y);

            if(ic<-marge || jc<-marge)
            continue;

            for(int qi=-1; qi<=2; ++qi)
            for(int qj=-p->j_dir; qj<=2*p->j_dir; ++qj)
            {
                const int ii=ic+qi, jj=jc+qj;

                if(ii<0 || ii>=p->knox || jj<0 || jj>=p->knoy)
                continue;

                // position along the ray and distance from it; for a sample slid along a wall, from the sample
                const double cx = p->XP[ii+marge] - (slid ? x : R.x);
                const double cy = (p->j_dir==1) ? p->YP[jj+marge] - (slid ? y : R.y) : 0.0;
                const double sc = (slid ? s : 0.0) + cx*R.nx + cy*R.ny;
                const double pd = fabs(-cx*R.ny + cy*R.nx);
                // half width of the painted strip: 1.5 cells, so that the strips of neighbouring crest points overlap
                // and reach the cells next to the side walls
                const double rr = 1.5*MAX(p->DXN[ii+marge], (p->j_dir==1) ? p->DYN[jj+marge] : 0.0);

                if(pd>rr || sc<sf0 || sc>sf1)
                continue;

                const int n = (ii-p->imin)*p->jmax + (jj-p->jmin);

                if(p->wet[n]==0)
                continue;

                brk[n] = 1;

                if(sc>=-R.sback && sc<=R.sa)
                {
                    const double v = R.nu*bc3_ramp((sc+R.sback)/MAX(R.sa+R.sback,1.0e-12));

                    if(v>vb(ii,jj))
                    {
                        vb(ii,jj) = v;

                        if(kw!=nullptr)
                        (*kw)(ii,jj) = R.kw;
                    }
                }
            }
        }
        }
    }
}

void breaking_crest3d::log(lexer *p, double t, const char *file)
{
    if(rank!=0)
    return;

    if(!logf.is_open())
    {
        logf.open(file);
        logf<<"# t id x y nx ny eta C us B state H h_c h_t b eps nu L dpn"<<endl;
    }

    for(const bc3_point &R : cand)
    logf<<t<<" "<<R.id<<" "<<R.x<<" "<<R.y<<" "<<R.nx<<" "<<R.ny<<" "<<R.eta<<" "<<R.C<<" "<<R.us<<" "<<R.B<<" "<<R.state
        <<" "<<R.H<<" "<<R.wl<<" "<<R.ht<<" "<<R.b<<" "<<R.eps<<" "<<R.nu<<" "<<R.L<<" "<<R.dpn<<"\n";
}
