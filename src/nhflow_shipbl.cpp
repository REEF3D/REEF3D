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
Author: Ahmet Soydan
--------------------------------------------------------------------*/

#include"nhflow_shipbl.h"
#include<Eigen/Dense>
#include<algorithm>
#include<cmath>
#include<cstdarg>
#include<cstdint>
#include<cstdio>
#include<map>

namespace
{
typedef nhflow_shipbl::P2 P2;
typedef nhflow_shipbl::tri3 tri3;
typedef nhflow_shipbl::curve curve;
typedef std::vector<double> vec;

const double PI = 3.14159265358979323846;

// numpy.interp: linear, constant beyond the ends
double interp(double x, const vec &xp, const vec &fp)
{
    const int n = int(xp.size());
    if(n==1 || x<=xp[0])
    return fp[0];
    if(x>=xp[n-1])
    return fp[n-1];
    const int j = int(std::upper_bound(xp.begin(),xp.end(),x) - xp.begin());
    return fp[j-1] + (x-xp[j-1])/(xp[j]-xp[j-1])*(fp[j]-fp[j-1]);
}

// closed section polygons (y, z) of the surface at x = x0
void slice_x(const std::vector<tri3> &tri, double x0, std::vector<curve> &loops)
{
    loops.clear();

    for(bool hit=true; hit;)
    {
        hit=false;
        for(const tri3 &t : tri)
        if(t[0]==x0 || t[3]==x0 || t[6]==x0)
        {
            hit=true;
            x0 += 1.0e-9;
            break;
        }
    }

    struct seg {P2 p[2]; int64_t k[2][2];};
    std::vector<seg> sg;

    for(const tri3 &t : tri)
    {
        const double lo = std::min(t[0],std::min(t[3],t[6]));
        const double hi = std::max(t[0],std::max(t[3],t[6]));
        if(!(lo<x0 && hi>x0))
        continue;

        seg q;
        int cnt=0;
        for(int e=0; e<3; ++e)
        {
            const double *pa = &t[3*e];
            const double *pb = &t[3*((e+1)%3)];
            if((pa[0]-x0)*(pb[0]-x0)<0.0)
            {
                const double *p1 = pa[0]>pb[0] ? pb : pa;      // same end points in the same order for both triangles of an edge
                const double *p2 = pa[0]>pb[0] ? pa : pb;
                const double w = (x0-p1[0])/(p2[0]-p1[0]);
                if(cnt<2)
                {
                    q.p[cnt].y = p1[1] + w*(p2[1]-p1[1]);
                    q.p[cnt].z = p1[2] + w*(p2[2]-p1[2]);
                }
                ++cnt;
            }
        }
        if(cnt==2)
        {
            for(int e=0; e<2; ++e)
            {
                q.k[e][0] = int64_t(std::nearbyint(q.p[e].y/1.0e-8));
                q.k[e][1] = int64_t(std::nearbyint(q.p[e].z/1.0e-8));
            }
            sg.push_back(q);
        }
    }

    typedef std::pair<int64_t,int64_t> key;
    std::map<key,std::vector<std::pair<int,int> > > adj;
    for(int n=0; n<int(sg.size()); ++n)
    for(int e=0; e<2; ++e)
    adj[key(sg[n].k[e][0],sg[n].k[e][1])].push_back(std::make_pair(n,e));

    std::vector<char> used(sg.size(),0);
    for(int s0=0; s0<int(sg.size()); ++s0)
    {
        if(used[s0])
        continue;
        used[s0]=1;
        curve loop;
        loop.push_back(sg[s0].p[0]);
        loop.push_back(sg[s0].p[1]);
        key k(sg[s0].k[1][0],sg[s0].k[1][1]);
        const key k0(sg[s0].k[0][0],sg[s0].k[0][1]);
        while(k!=k0)
        {
            int sn=-1, en=0;
            for(const std::pair<int,int> &c : adj[k])
            if(!used[c.first])
            {
                sn=c.first;
                en=c.second;
                break;
            }
            if(sn<0)
            break;
            used[sn]=1;
            loop.push_back(sg[sn].p[1-en]);
            k = key(sg[sn].k[1-en][0],sg[sn].k[1-en][1]);
        }
        if(k==k0 && loop.size()>3)
        {
            loop.pop_back();
            loops.push_back(loop);
        }
    }
}

// the parts of a closed section polygon on the port side (y >= 0): open polylines from and to the centre plane
void port_parts(const curve &loop_in, std::vector<curve> &out, double eps=1.0e-9)
{
    const int n = int(loop_in.size());
    std::vector<int> sgn(n);
    int first=-1;
    for(int k=0; k<n; ++k)
    {
        sgn[k] = loop_in[k].y>eps ? 1 : (loop_in[k].y<-eps ? -1 : 0);
        if(sgn[k]<0 && first<0)
        first=k;
    }
    if(first<0)
    return;

    auto cross = [](const P2 &a, const P2 &b)
    {
        const double w = a.y/(a.y-b.y);
        P2 r;
        r.y = 0.0;
        r.z = a.z + w*(b.z-a.z);
        return r;
    };

    curve cur;
    bool open=false;
    for(int kk=0; kk<n; ++kk)
    {
        const int k = (kk+first)%n;
        const int k1 = (kk+1+first)%n;
        const P2 &a = loop_in[k];
        const P2 &b = loop_in[k1];
        const int sa = sgn[k], sb = sgn[k1];

        if(sa<0)
        {
            if(sb>=0)
            {
                cur.clear();
                open=true;
                if(sb!=0)
                cur.push_back(cross(a,b));
            }
            continue;
        }
        if(!open)
        {
            cur.clear();
            open=true;
        }
        P2 q;
        q.y = std::max(a.y,0.0);
        q.z = a.z;
        cur.push_back(q);
        if(sb<0)
        {
            if(sa>0)
            cur.push_back(cross(a,b));
            curve p = cur;
            cur.clear();
            open=false;
            double ymax=-1.0e30;
            for(const P2 &v : p)
            ymax = std::max(ymax,v.y);
            if(p.size()>1 && ymax>1.0e-7)
            {
                if(p.front().z>p.back().z)
                std::reverse(p.begin(),p.end());
                p.front().y = 0.0;
                p.back().y = 0.0;
                out.push_back(p);
            }
        }
    }
}

// one curve from (0, z_env) through all hull parts to (0, z_cap): points and the wall flag of the segments
void station_curve(std::vector<curve> parts, double z_env, double z_cap, curve &pts, std::vector<char> &wl)
{
    std::stable_sort(parts.begin(),parts.end(),[](const curve &a, const curve &b){return a[0].z<b[0].z;});
    pts.clear();
    wl.clear();
    P2 q;
    q.y=0.0;
    q.z=z_env;
    pts.push_back(q);
    for(const curve &p : parts)
    {
        if(fabs(p[0].z-pts.back().z)>1.0e-9 || fabs(pts.back().y)>1.0e-12)
        {
            pts.push_back(p[0]);
            wl.push_back(0);
        }
        for(size_t m=1; m<p.size(); ++m)
        {
            pts.push_back(p[m]);
            wl.push_back(1);
        }
    }
    if(z_cap-pts.back().z>1.0e-9)
    {
        q.z=z_cap;
        pts.push_back(q);
        wl.push_back(0);
    }
}

// the part of the curve up to its first crossing of the still water level z = 0
void cut_at_water(curve &pts, std::vector<char> &wl)
{
    for(size_t k=0; k+1<pts.size(); ++k)
    if(pts[k].z<0.0 && 0.0<=pts[k+1].z)
    {
        const double w = -pts[k].z/(pts[k+1].z-pts[k].z);
        P2 q;
        q.y = pts[k].y + w*(pts[k+1].y-pts[k].y);
        q.z = pts[k].z + w*(pts[k+1].z-pts[k].z);
        pts.resize(k+1);
        pts.push_back(q);
        wl.resize(k+1);
        return;
    }
}

void arcs_of(const curve &pts, vec &a)
{
    a.assign(pts.size(),0.0);
    for(size_t k=1; k<pts.size(); ++k)
    a[k] = a[k-1] + sqrt(pow(pts[k].y-pts[k-1].y,2.0) + pow(pts[k].z-pts[k-1].z,2.0));
}

// points of the curve at the fractions frac of its arc length, the wall flag there, and the arc length
double resample(const curve &pts, const std::vector<char> &wl, const vec &frac, P2 *q, char *w)
{
    vec a, py(pts.size()), pz(pts.size());
    arcs_of(pts,a);
    for(size_t k=0; k<pts.size(); ++k)
    {
        py[k]=pts[k].y;
        pz[k]=pts[k].z;
    }
    const int nd = int(pts.size())-1;
    for(size_t m=0; m<frac.size(); ++m)
    {
        if(a.back()<=0.0)
        {
            q[m]=pts[0];
            if(w!=nullptr)
            w[m]=0;
            continue;
        }
        const double v = frac[m]*a.back();
        q[m].y = interp(v,a,py);
        q[m].z = interp(v,a,pz);
        if(w!=nullptr)
        {
            int k = int(std::upper_bound(a.begin(),a.end(),v) - a.begin()) - 1;
            k = std::max(0,std::min(k,nd-1));
            w[m] = wl[k];
        }
    }
    return std::max(a.back(),0.0);
}

// points of the curve at the arc lengths arcs from its start (beyond the end: the end point)
double resample_abs(const curve &pts, const vec &arcs, curve &q)
{
    vec a, py(pts.size()), pz(pts.size());
    arcs_of(pts,a);
    for(size_t k=0; k<pts.size(); ++k)
    {
        py[k]=pts[k].y;
        pz[k]=pts[k].z;
    }
    q.resize(arcs.size());
    for(size_t m=0; m<arcs.size(); ++m)
    {
        const double v = std::min(arcs[m],a.back());
        q[m].y = interp(v,a,py);
        q[m].z = interp(v,a,pz);
    }
    return a.back();
}

double H1_of(double H)
{
    H = std::max(1.12,std::min(H,3.0));
    return H<=1.6 ? 3.3 + 0.8234*pow(H-1.1,-1.287) : 3.3 + 1.5501*pow(H-0.6778,-3.064);
}

double H_of(double H1)
{
    H1 = std::max(H1,3.32);
    return H1>=5.3 ? 1.1 + pow((H1-3.3)/0.8234,-1.0/1.287) : 0.6778 + pow((H1-3.3)/1.5501,-1.0/3.064);
}

double taper(double s, double s_te)
{
    const double r = s-s_te;
    return r<0.5 ? 1.0 : (r<0.9 ? 0.5*(1.0+cos(PI*(r-0.5)/0.4)) : 0.0);
}

// three-point smoothing along the stations, n passes; f is (ns, nz)
void smooth(vec &f, int ns, int nz, int n)
{
    vec g(f.size());
    for(int q=0; q<n; ++q)
    {
        g=f;
        for(int i=1; i<ns-1; ++i)
        for(int k=0; k<nz; ++k)
        g[i*nz+k] = 0.25*f[(i-1)*nz+k] + 0.5*f[i*nz+k] + 0.25*f[(i+1)*nz+k];
        f=g;
    }
}
}

nhflow_shipbl::nhflow_shipbl() : L(1.0),T(0.0),B(0.0),U(1.0),nu(1.0e-6),x_fp(0.0),y_c(0.0),z_w(0.0),
                                 xtr(0.05),cwake(0.3),cdiff(0.04),fsw(true),ds(0.00625),s_end(1.9),nz(24),
                                 CF(0.0),CV(0.0),Sw(0.0),xmin(0.0),xmax(0.0),ns(0)
{
}

void nhflow_shipbl::say(const char *fmt, ...)
{
    char buf[1024];
    va_list ap;
    va_start(ap,fmt);
    vsnprintf(buf,sizeof(buf),fmt,ap);
    va_end(ap);
    log += buf;
    log += "\n";
}

void nhflow_shipbl::dump(const char *name, const std::vector<double> &f) const
{
    if(dumpdir.empty())
    return;
    FILE *fp = fopen((dumpdir+"/"+name+".bin").c_str(),"wb");
    if(fp==nullptr)
    return;
    fwrite(f.data(),sizeof(double),f.size(),fp);
    fclose(fp);
}

void nhflow_shipbl::dump(const char *name, const std::vector<P2> &f) const
{
    vec g(2*f.size());
    for(size_t n=0; n<f.size(); ++n)
    {
        g[2*n]=f[n].y;
        g[2*n+1]=f[n].z;
    }
    dump(name,g);
}

// ---------------------------------------------------------------- sections of the hull
void nhflow_shipbl::set_hull(const std::vector<tri3> &hull)
{
    tri = hull;
    xmin=1.0e30;
    xmax=-1.0e30;
    double ylo=1.0e30, yhi=-1.0e30;
    for(tri3 &t : tri)
    for(int v=0; v<3; ++v)
    {
        t[3*v] -= x_fp;
        t[3*v+1] -= y_c;
        t[3*v+2] -= z_w;
        xmin = std::min(xmin,t[3*v]);
        xmax = std::max(xmax,t[3*v]);
        ylo = std::min(ylo,t[3*v+1]);
        yhi = std::max(yhi,t[3*v+1]);
    }
    if(B<=0.0)
    B = yhi-ylo;

    // depth envelope: deepest point of the hull up to x
    env_x.resize(799);
    env_z.resize(799);
    std::vector<curve> loops, pp;
    for(int m=0; m<799; ++m)
    {
        env_x[m] = xmin + (xmax-xmin)*double(m+1)/800.0;
        slice_x(tri,env_x[m],loops);
        pp.clear();
        for(const curve &l : loops)
        port_parts(l,pp);
        double zb = pp.empty() ? 0.0 : 1.0e30;
        for(const curve &p : pp)
        for(const P2 &v : p)
        zb = std::min(zb,v.z);
        env_z[m] = m>0 ? std::min(zb,env_z[m-1]) : zb;
    }
    if(T<=0.0)
    T = -env_z.back();
}

double nhflow_shipbl::z_env(double x) const
{
    return interp(x,env_x,env_z);
}

void nhflow_shipbl::parts(double x, std::vector<curve> &pp) const
{
    pp.clear();
    if(x<=xmin || x>=xmax)
    return;
    std::vector<curve> loops;
    slice_x(tri,x,loops);
    for(const curve &l : loops)
    port_parts(l,pp);
}

void nhflow_shipbl::wet(double x, curve &pts, std::vector<char> &wl) const
{
    std::vector<curve> pp;
    parts(x,pp);
    station_curve(pp,z_env(x),0.0,pts,wl);
    cut_at_water(pts,wl);
}

// to the deck (hull stations) or to z_top (stations without hull above the water)
void nhflow_shipbl::full(double x, double z_top, curve &pts, std::vector<char> &wl) const
{
    std::vector<curve> pp;
    parts(x,pp);
    double ztop = z_top;
    for(const curve &p : pp)
    ztop = std::max(ztop,p.back().z);
    station_curve(pp,z_env(x),ztop,pts,wl);
}

// ---------------------------------------------------------------- surface grid
void nhflow_shipbl::build_grid()
{
    const int k0 = int(floor(xmin/L/ds)) - 1;
    const int k1 = int(std::nearbyint(s_end/ds));
    ns = k1-k0+1;
    s.resize(ns);
    for(int i=0; i<ns; ++i)
    s[i] = double(k0+i)*ds;

    vec fn(nz+1), fc(nz);
    for(int k=0; k<=nz; ++k)
    fn[k] = double(k)/double(nz);
    for(int k=0; k<nz; ++k)
    fc[k] = (double(k)+0.5)/double(nz);

    Pn.resize(ns*(nz+1));
    Pc.resize(ns*nz);
    wall.resize(ns*nz);
    girth.resize(ns);
    Hn.resize((ns+1)*(nz+1));
    Hc.resize((ns+1)*nz);

    curve p;
    std::vector<char> w;
    for(int i=0; i<ns; ++i)
    {
        wet(s[i]*L,p,w);
        resample(p,w,fn,&Pn[i*(nz+1)],nullptr);
        girth[i] = resample(p,w,fc,&Pc[i*nz],&wall[i*nz]);
    }
    for(int i=0; i<=ns; ++i)
    {
        wet((s[0] + (double(i)-0.5)*ds)*L,p,w);
        resample(p,w,fn,&Hn[i*(nz+1)],nullptr);
        resample(p,w,fc,&Hc[i*nz],nullptr);
    }
    metrics();

    int nw=0;
    double A=0.0, gmax=0.0;
    for(int i=0; i<ns; ++i)
    {
        gmax = std::max(gmax,girth[i]);
        for(int k=0; k<nz; ++k)
        if(wall[i*nz+k])
        {
            ++nw;
            A += J[i*nz+k];
        }
    }
    A *= 2.0*ds*L/double(nz);
    say("ship boundary layer: %d stations s = %.4f .. %.3f (hull %.4f .. %.4f), %d cells along the girth, wall cells %d",
        ns,s[0],s[ns-1],xmin/L,xmax/L,nz,nw);
    say("  wetted surface of the grid %.4f m^2; largest girth %.3f m, draft %.4f m",A,gmax,-env_z.back());

    dump("Pn",Pn);
    dump("Pc",Pc);
    dump("Hn",Hn);
    dump("Hc",Hc);
    dump("girth",girth);
    dump("wall",vec(wall.begin(),wall.end()));
    dump("J",J);
}

void nhflow_shipbl::metrics()
{
    const double dx = ds*L;
    const int N = ns*nz;
    e1y.resize(N); e1z.resize(N); e2y.resize(N); e2z.resize(N);
    nx.resize(N); ny.resize(N); nzc.resize(N); J.resize(N); h2.resize(N);

    for(int i=0; i<ns; ++i)
    for(int k=0; k<nz; ++k)
    {
        const int m = i*nz+k;
        const double yx = std::max(-5.0,std::min((Hc[(i+1)*nz+k].y - Hc[i*nz+k].y)/dx,5.0));
        const double zx = std::max(-5.0,std::min((Hc[(i+1)*nz+k].z - Hc[i*nz+k].z)/dx,5.0));
        const double yt = (Pn[i*(nz+1)+k+1].y - Pn[i*(nz+1)+k].y)*double(nz);
        double zt = (Pn[i*(nz+1)+k+1].z - Pn[i*(nz+1)+k].z)*double(nz);
        if(yt*yt + zt*zt < 1.0e-16)
        zt = 1.0e-6;
        e1y[m]=yx; e1z[m]=zx; e2y[m]=yt; e2z[m]=zt;
        // N = e2 x e1 with e1 = (1, yx, zx), e2 = (0, yt, zt)
        const double Nx = yt*zx - zt*yx, Ny = zt, Nz = -yt;
        J[m] = sqrt(Nx*Nx + Ny*Ny + Nz*Nz);
        nx[m]=Nx/J[m]; ny[m]=Ny/J[m]; nzc[m]=Nz/J[m];
        h2[m] = sqrt(yt*yt + zt*zt);
    }
}

// ---------------------------------------------------------------- double body (Hess & Smith, port side unknowns)
void nhflow_shipbl::double_body()
{
    struct V3 {double x,y,z;};
    std::vector<std::array<V3,3> > pt;
    std::vector<int> own;

    for(int i=0; i<ns; ++i)
    for(int k=0; k<nz; ++k)
    {
        if(!wall[i*nz+k])
        continue;
        const double xa = (s[0] + (double(i)-0.5)*ds)*L, xb = (s[0] + (double(i)+0.5)*ds)*L;
        const V3 a = {xa,Hn[i*(nz+1)+k].y,Hn[i*(nz+1)+k].z};
        const V3 b = {xb,Hn[(i+1)*(nz+1)+k].y,Hn[(i+1)*(nz+1)+k].z};
        const V3 c = {xa,Hn[i*(nz+1)+k+1].y,Hn[i*(nz+1)+k+1].z};
        const V3 d = {xb,Hn[(i+1)*(nz+1)+k+1].y,Hn[(i+1)*(nz+1)+k+1].z};
        const std::array<V3,3> two[2] = {{{a,d,b}},{{a,c,d}}};
        for(int q=0; q<2; ++q)
        {
            const V3 &p0=two[q][0], &p1=two[q][1], &p2=two[q][2];
            const double cx = (p1.y-p0.y)*(p2.z-p0.z) - (p1.z-p0.z)*(p2.y-p0.y);
            const double cy = (p1.z-p0.z)*(p2.x-p0.x) - (p1.x-p0.x)*(p2.z-p0.z);
            const double cz = (p1.x-p0.x)*(p2.y-p0.y) - (p1.y-p0.y)*(p2.x-p0.x);
            if(sqrt(cx*cx + cy*cy + cz*cz)>1.0e-12)
            {
                pt.push_back(two[q]);
                own.push_back(i*nz+k);
            }
        }
    }

    const int N = int(pt.size());
    vec A(N), pnx(N), pny(N), pnz(N), cen(3*N), sub(12*N);
    for(int n=0; n<N; ++n)
    {
        const V3 &p0=pt[n][0], &p1=pt[n][1], &p2=pt[n][2];
        const double cx = (p1.y-p0.y)*(p2.z-p0.z) - (p1.z-p0.z)*(p2.y-p0.y);
        const double cy = (p1.z-p0.z)*(p2.x-p0.x) - (p1.x-p0.x)*(p2.z-p0.z);
        const double cz = (p1.x-p0.x)*(p2.y-p0.y) - (p1.y-p0.y)*(p2.x-p0.x);
        const double l = sqrt(cx*cx + cy*cy + cz*cz);
        A[n] = 0.5*l;
        pnx[n]=cx/l; pny[n]=cy/l; pnz[n]=cz/l;
        const double v[3][3] = {{p0.x,p0.y,p0.z},{p1.x,p1.y,p1.z},{p2.x,p2.y,p2.z}};
        for(int q=0; q<3; ++q)
        {
            cen[3*n+q] = (v[0][q]+v[1][q]+v[2][q])/3.0;
            const double m0 = 0.5*(v[0][q]+v[1][q]), m1 = 0.5*(v[1][q]+v[2][q]), m2 = 0.5*(v[2][q]+v[0][q]);
            sub[12*n+q]   = (v[0][q]+m0+m2)/3.0;
            sub[12*n+3+q] = (m0+v[1][q]+m1)/3.0;
            sub[12*n+6+q] = (m2+m1+v[2][q])/3.0;
            sub[12*n+9+q] = (m0+m1+m2)/3.0;
        }
    }
    const double img[4][3] = {{1.0,1.0,1.0},{1.0,-1.0,1.0},{1.0,1.0,-1.0},{1.0,-1.0,-1.0}};

    Eigen::MatrixXd Am(N,N);
    Eigen::VectorXd rhs(N);
    for(int i=0; i<N; ++i)
    {
        rhs(i) = -pnx[i];
        for(int j=0; j<N; ++j)
        {
            double acc=0.0;
            for(int im=0; im<4; ++im)
            {
                if(im==0 && i==j)
                continue;
                for(int q=0; q<4; ++q)
                {
                    const double dx = cen[3*i]   - sub[12*j+3*q]*img[im][0];
                    const double dy = cen[3*i+1] - sub[12*j+3*q+1]*img[im][1];
                    const double dz = cen[3*i+2] - sub[12*j+3*q+2]*img[im][2];
                    const double r2 = std::max(dx*dx + dy*dy + dz*dz,1.0e-30);
                    acc += (dx*pnx[i] + dy*pny[i] + dz*pnz[i])/(4.0*PI*r2*sqrt(r2));
                }
            }
            Am(i,j) = acc*0.25*A[j];
        }
        Am(i,i) += 0.5;
    }
    const Eigen::VectorXd sig = Am.partialPivLu().solve(rhs);
    Am.resize(0,0);

    // perturbation velocity at a point of the sources on the port panels and their three images
    auto vel = [&](double px, double py, double pz, int skip, double *out)
    {
        double ox=0.0, oy=0.0, oz=0.0;
        for(int j=0; j<N; ++j)
        {
            const double qj = 0.25*sig(j)*A[j]/(4.0*PI);
            for(int im=0; im<4; ++im)
            {
                if(im==0 && j==skip)
                continue;
                for(int q=0; q<4; ++q)
                {
                    const double dx = px - sub[12*j+3*q]*img[im][0];
                    const double dy = py - sub[12*j+3*q+1]*img[im][1];
                    const double dz = pz - sub[12*j+3*q+2]*img[im][2];
                    const double r2 = std::max(dx*dx + dy*dy + dz*dz,1.0e-30);
                    const double w = qj/(r2*sqrt(r2));
                    ox += dx*w; oy += dy*w; oz += dz*w;
                }
            }
        }
        out[0]=ox; out[1]=oy; out[2]=oz;
    };

    const int M = ns*nz;
    Vx.assign(M,0.0); Vy.assign(M,0.0); Vz.assign(M,0.0);
    std::vector<int> cnt(M,0);
    double vnmax=0.0, ssum=0.0, asum=0.0;
    for(int n=0; n<N; ++n)
    {
        double v[3];
        vel(cen[3*n],cen[3*n+1],cen[3*n+2],n,v);
        v[0] += 0.5*sig(n)*pnx[n] + 1.0;
        v[1] += 0.5*sig(n)*pny[n];
        v[2] += 0.5*sig(n)*pnz[n];
        vnmax = std::max(vnmax,fabs(v[0]*pnx[n] + v[1]*pny[n] + v[2]*pnz[n]));
        ssum += sig(n)*A[n];
        asum += A[n];
        Vx[own[n]] += v[0]; Vy[own[n]] += v[1]; Vz[own[n]] += v[2];
        ++cnt[own[n]];
    }
    for(int i=0; i<ns; ++i)
    for(int k=0; k<nz; ++k)
    {
        const int m = i*nz+k;
        if(cnt[m]>0)
        {
            Vx[m] /= double(cnt[m]); Vy[m] /= double(cnt[m]); Vz[m] /= double(cnt[m]);
        }
        else                                                   // no wall: the point on the centre plane itself
        {
            double v[3];
            vel(s[i]*L,std::max(Pc[m].y,0.0),std::min(Pc[m].z,-1.0e-6),-1,v);
            Vx[m]=v[0]+1.0; Vy[m]=v[1]; Vz[m]=v[2];
        }
    }
    say("  double body: %d panels (port side), wetted area %.4f m^2, sum of the sources %.2e, largest normal velocity at the centroids %.1e",
        N,2.0*asum,ssum,vnmax);
    dump("Vx",Vx);
    dump("Vy",Vy);
    dump("Vz",Vz);
}

// ---------------------------------------------------------------- integral boundary layer
// Conservation form in (x, t): with W = J a/Ue (width of the stream tube of a cell per unit t) and q2 = c/Ue
//     d(W theta)/dx = J [cf/2 - (H + 2) theta/Ue dUe/ds] - d(J q2 theta)/dt  (+ girthwise diffusion)
//     d(W E)/dx     = J Ue F                             - d(J q2 E)/dt,      E = Ue theta H1
// No flux through t = 0: the centre line of a flat bottom (symmetry) or a sharp keel.
void nhflow_shipbl::solve()
{
    const int M = ns*nz;
    const double dx = ds*L, dt_ = 1.0/double(nz);
    const int nsm = 8;

    vec t0(M), t1(M), t2(M), a(M), c(M);
    for(int m=0; m<M; ++m)
    {
        const double vx=U*Vx[m], vy=U*Vy[m], vz=U*Vz[m];
        const double vn = vx*nx[m] + vy*ny[m] + vz*nzc[m];
        t0[m] = vx - vn*nx[m];
        t1[m] = vy - vn*ny[m];
        t2[m] = vz - vn*nzc[m];
        // c = (Vt . e2 - Vt_x e1 . e2)/h2^2
        c[m] = ((t1[m]*e2y[m] + t2[m]*e2z[m]) - t0[m]*(e1y[m]*e2y[m] + e1z[m]*e2z[m]))/(h2[m]*h2[m]);
    }
    smooth(c,ns,nz,nsm);
    smooth(t0,ns,nz,nsm);
    smooth(t1,ns,nz,nsm);
    smooth(t2,ns,nz,nsm);
    a = t0;

    Ue.resize(M);
    vec am(M), W(M), Q2(M), dUds(M);
    for(int m=0; m<M; ++m)
    {
        Ue[m] = std::max(sqrt(t0[m]*t0[m] + t1[m]*t1[m] + t2[m]*t2[m]),0.05*U);
        am[m] = std::max(a[m],0.3*U);                          // marching speed (stagnation regions)
        W[m] = J[m]*am[m]/Ue[m];
        Q2[m] = J[m]*c[m]/Ue[m];
    }
    for(int i=0; i<ns; ++i)
    for(int k=0; k<nz; ++k)
    {
        const int m = i*nz+k;
        const double Uex = i==0 ? (Ue[m+nz]-Ue[m])/dx : (i==ns-1 ? (Ue[m]-Ue[m-nz])/dx : (Ue[m+nz]-Ue[m-nz])/(2.0*dx));
        const double Uet = k==0 ? (Ue[m+1]-Ue[m])/dt_ : (k==nz-1 ? (Ue[m]-Ue[m-1])/dt_ : (Ue[m+1]-Ue[m-1])/(2.0*dt_));
        dUds[m] = (am[m]*Uex + c[m]*Uet)/Ue[m];
    }

    th.assign(M,1.0e-6);
    H.assign(M,1.4);
    cf.assign(M,0.0);
    Ms.assign(ns,0.0);
    Ds.assign(ns,0.0);

    // laminar part (Thwaites) along the grid lines, from where the wall starts to the trip
    const int itr = int(std::nearbyint((xtr-s[0])/ds));
    vec I(nz,0.0);
    for(int i=0; i<=itr; ++i)
    for(int k=0; k<nz; ++k)
    {
        const int m = i*nz+k;
        if(wall[m])
        {
            I[k] += pow(Ue[m],5.0)*(Ue[m]/am[m])*dx;
            th[m] = std::max(sqrt(0.45*nu*I[k]/pow(Ue[m],6.0)),1.0e-6);
            H[m] = 2.59;
            cf[m] = 2.0*0.22*nu/(Ue[m]*std::max(th[m],1.0e-6));
        }
        else
        I[k] = 0.0;
    }

    vec t(nz), E(nz), G(nz,pow(1.4-1.0,-2.0));
    for(int k=0; k<nz; ++k)
    {
        t[k] = th[itr*nz+k];
        E[k] = Ue[itr*nz+k]*t[k]*H1_of(1.4);
    }

    // upwind flux J q2 f through the faces between the cells (nz + 1 values)
    auto flux = [&](const vec &q2, const vec &f, bool w0, vec &out)
    {
        out[0] = (w0 ? 0.0 : 1.5*q2[0]-0.5*q2[1])*f[0];
        for(int k=1; k<nz; ++k)
        {
            const double qf = 0.5*(q2[k]+q2[k-1]);
            out[k] = qf>0.0 ? qf*f[k-1] : qf*f[k];
        }
        out[nz] = (1.5*q2[nz-1]-0.5*q2[nz-2])*f[nz-1];
    };
    auto dif = [&](const vec &f, const vec &dl, const vec &w, const vec &hh, vec &out)
    {
        vec fl(nz-1);
        for(int k=0; k<nz-1; ++k)
        fl[k] = 0.5*(w[k+1]*dl[k+1]/(hh[k+1]*hh[k+1]) + w[k]*dl[k]/(hh[k]*hh[k]))*(f[k+1]-f[k])/dt_;
        out[0] = cdiff*fl[0]/dt_;
        for(int k=1; k<nz-1; ++k)
        out[k] = cdiff*(fl[k]-fl[k-1])/dt_;
        out[nz-1] = -cdiff*fl[nz-2]/dt_;
    };

    double ms=0.0, hm=0.0, sth=0.0, sab=0.0;
    int nsubmax=0;
    vec hh2(nz), m_(nz), ME(nz), Wf(nz), Qf(nz), UE(nz), Jf(nz), Af(nz), Cf(nz), dU(nz), hw(nz), H1(nz), hk(nz), hc(nz),
        cfl(nz), F(nz), dl(nz), ft(nz+1), fE(nz+1), d1(nz), d2(nz), d3(nz), d4(nz), hct(nz), Wn(nz), UEn(nz);

    for(int i=itr; i<ns-1; ++i)
    {
        const char *w = &wall[(i+1)*nz];
        double qmax=0.0, dmax=0.0;
        for(int k=0; k<nz; ++k)
        {
            const int m = i*nz+k;
            hh2[k] = 0.5*(h2[m]+h2[m+nz]);
            const double dl0 = t[k]*std::min(H1_of(std::max(1.12,std::min(H[m],3.0))) + H[m],15.0);
            qmax = std::max(qmax,fabs(0.5*(Q2[m]/W[m] + Q2[m+nz]/W[m+nz])));
            dmax = std::max(dmax,dl0/pow(hh2[k]*dt_,2.0));
        }
        const int nsub = int(std::min(std::max(5.0,std::max(ceil(1.3*qmax*dx/dt_),ceil(3.0*cdiff*dx*dmax))),400.0));
        nsubmax = std::max(nsubmax,nsub);
        const double h_ = dx/double(nsub);
        for(int k=0; k<nz; ++k)
        {
            m_[k] = W[i*nz+k]*t[k];
            ME[k] = W[i*nz+k]*E[k];
        }

        for(int q=0; q<nsub; ++q)
        {
            const double f = (double(q)+0.5)/double(nsub), f1 = (double(q)+1.0)/double(nsub);
            for(int k=0; k<nz; ++k)
            {
                const int m = i*nz+k;
                Wf[k] = (1.0-f)*W[m] + f*W[m+nz];
                Qf[k] = (1.0-f)*Q2[m] + f*Q2[m+nz];
                UE[k] = (1.0-f)*Ue[m] + f*Ue[m+nz];
                Jf[k] = (1.0-f)*J[m] + f*J[m+nz];
                Af[k] = (1.0-f)*am[m] + f*am[m+nz];
                Cf[k] = (1.0-f)*c[m] + f*c[m+nz];
                dU[k] = (1.0-f)*dUds[m] + f*dUds[m+nz];
                Wn[k] = (1.0-f1)*W[m] + f1*W[m+nz];
                UEn[k] = (1.0-f1)*Ue[m] + f1*Ue[m+nz];
                t[k] = m_[k]/Wf[k];
                E[k] = ME[k]/Wf[k];
                H1[k] = E[k]/(UE[k]*t[k]);
                hw[k] = H_of(H1[k]);
                hk[k] = 1.0 + pow(G[k],-0.5);
                hc[k] = w[k] ? hw[k] : hk[k];
                cfl[k] = w[k] ? 0.246*pow(10.0,-0.678*hw[k])*pow(std::max(UE[k]*t[k]/nu,50.0),-0.268) : 0.0;
                F[k] = 0.0306*pow(std::max(H1[k]-3.0,0.3),-0.6169);
                dl[k] = t[k]*std::min((w[k] ? H1[k] : H1_of(hk[k])) + hc[k],15.0);       // boundary-layer thickness
                hct[k] = hc[k]*t[k];
            }
            flux(Qf,t,w[0]!=0,ft);
            flux(Qf,E,w[0]!=0,fE);
            ms += h_*ft[nz];
            hm += h_*ft[nz]*hc[nz-1];
            sth += fabs(h_*ft[nz])*t[nz-1];
            sab += fabs(h_*ft[nz]);
            dif(t,dl,Wf,hh2,d1);
            dif(E,dl,Wf,hh2,d2);
            dif(hct,dl,Wf,hh2,d3);
            for(int k=0; k<nz; ++k)
            {
                const double dm = Jf[k]*(0.5*cfl[k] - (hc[k]+2.0)*t[k]/UE[k]*dU[k]) - (ft[k+1]-ft[k])/dt_ + d1[k];
                const double dME = Jf[k]*UE[k]*F[k] - (fE[k+1]-fE[k])/dt_ + d2[k];
                const double lo = k==0 ? 0.0 : (G[k]-G[k-1])/dt_;
                const double hi = k==nz-1 ? 0.0 : (G[k+1]-G[k])/dt_;
                d4[k] = (UE[k]*cwake/t[k] - Cf[k]*(Cf[k]>0.0 ? lo : hi))/Af[k];
                const double dsn = (Wf[k]*hc[k]*t[k] + h_*d3[k])/Wn[k];          // delta* after the diffusion (wake)
                m_[k] = std::max(m_[k] + h_*dm,1.0e-6*Wn[k]);
                ME[k] = std::max(ME[k] + h_*dME,3.32*UEn[k]*m_[k]);
                d3[k] = dsn;
            }
            for(int k=0; k<nz; ++k)
            {
                double Gn = G[k] + h_*d4[k];
                if(cdiff>0.0)
                {
                    const double hn = std::max(1.005,std::min(d3[k]/(Wf[k]*t[k]/Wn[k] + h_*d1[k]/Wn[k]),3.0));
                    Gn = std::max(Gn + (pow(hn-1.0,-2.0) - pow(hk[k]-1.0,-2.0)),0.1);
                }
                // keep the variable of the other kind consistent (a cell can change its kind from station to station)
                const double hw_n = H_of(ME[k]/(UEn[k]*m_[k]));
                const double hk_n = 1.0 + pow(Gn,-0.5);
                if(!w[k])
                ME[k] = UEn[k]*m_[k]*H1_of(std::max(hk_n,1.12));
                d4[k] = w[k] ? pow(std::max(hw_n,1.02)-1.0,-2.0) : Gn;
            }
            G = d4;
        }

        for(int k=0; k<nz; ++k)
        {
            const int m = (i+1)*nz+k;
            t[k] = m_[k]/W[m];
            E[k] = ME[k]/W[m];
            const double hh = w[k] ? H_of(E[k]/(Ue[m]*t[k])) : 1.0 + pow(G[k],-0.5);
            th[m] = t[k];
            H[m] = hh;
            cf[m] = w[k] ? 0.246*pow(10.0,-0.678*hh)*pow(std::max(Ue[m]*t[k]/nu,50.0),-0.268) : 0.0;
        }
        if(ms>1.0e-9)
        {
            double hs = std::min(std::max(hm/ms,1.02),3.0);
            ms *= pow(Ue[i*nz+nz-1]/Ue[(i+1)*nz+nz-1],hs+2.0);     // momentum equation of a wake in the pressure gradient at the surface
            hs = 1.0 + pow(pow(hs-1.0,-2.0) + cwake*dx/std::max(sth/std::max(sab,1.0e-30),1.0e-4),-0.5);
            hm = hs*ms;
            Ms[i+1] = ms;
            Ds[i+1] = hm;
        }
    }

    // friction and viscous resistance of the layer (Squire-Young in the wake + surface wake)
    const double Re = U*L/nu, te = xmax/L;
    double tx=0.0;
    Sw=0.0;
    for(int m=0; m<M; ++m)
    if(wall[m])
    {
        const double dA = J[m]*dx*dt_;
        Sw += 2.0*dA;
        tx += 2.0*cf[m]*pow(Ue[m]/U,2.0)*(std::max(a[m],0.0)/Ue[m])*dA;
    }
    CF = tx/Sw;
    const int j0 = int(std::nearbyint((te+0.3-s[0])/ds)), j1 = int(std::nearbyint((te+0.5-s[0])/ds));
    const int dte = int(std::nearbyint((te-s[0])/ds));
    double thi=0.0, msw=0.0;
    for(int i=j0; i<j1; ++i)
    {
        for(int k=0; k<nz; ++k)
        {
            const int m = i*nz+k;
            thi += th[m]*pow(Ue[m]/U,(H[m]+5.0)/2.0)*W[m]*dt_;
        }
        msw += Ms[i]*pow(Ue[i*nz+nz-1]/U,3.1);
    }
    thi /= double(j1-j0);
    msw /= double(j1-j0);
    CV = 4.0*(thi+msw)/Sw;
    const double cfi = 0.075/pow(log10(Re)-2.0,2.0);

    dstar.resize(M);
    dplane.resize(M);
    for(int i=0; i<ns; ++i)
    {
        const double tp = taper(s[i],te);
        for(int k=0; k<nz; ++k)
        {
            const int m = i*nz+k;
            dstar[m] = H[m]*th[m];
            dplane[m] = dstar[m]*std::min(J[m]/h2[m],3.0)*tp;          // displacement in the section plane
        }
        Ds[i] = std::max(Ds[i],0.0)*tp;
    }
    say("  boundary layer: Re %.2e, U %.4f m/s, trip at s = %.3f, cross transport up to %d substeps",Re,U,xtr,nsubmax);
    say("  surface wake (left through the still water line): momentum area %.1f cm^2 of %.1f in all (one side), displacement area at the stern end %.1f cm^2",
        1.0e4*msw,1.0e4*(thi+msw),1.0e4*Ds[dte]);
    say("  wetted surface of the grid %.4f m^2; C_F = %.3f e-3 (ITTC-57 %.3f); C_V from the wake %.3f e-3 -> C_V / C_F = %.3f, C_V / C_F,ITTC = %.3f",
        Sw,1.0e3*CF,1.0e3*cfi,1.0e3*CV,CV/CF,CV/cfi);
    dump("theta",th);
    dump("H",H);
    dump("dstar",dstar);
    dump("dplane",dplane);
    dump("Ds",Ds);
    dump("Ms",Ms);
}

// ---------------------------------------------------------------- displacement body
// Structured surface, port and starboard. Nodes at fixed arc lengths from the start of every section curve. Above the
// still water level: the displacement of the top row over the arc w from the water line, then down to 0 at 1.5 w;
// w = 0.3 T n_y^2 (steep side: the layer under a wave crest) + Ds/(1.25 delta*) (flat overhang: the surface wake Ds
// that left through the water line lies under the dry part of the overhang). Behind the hull the surface wake is a
// thickening of the centre-plane wake, largest at the surface and 0 at the depth T (area Ds).
void nhflow_shipbl::body(std::vector<tri3> &out)
{
    const int nst=300, ntail=100, nn=110;
    const double fb = 0.5*T, eps = 1.0e-4*L;

    vec tc(nz), row0(nz), row1(nz), Dsw(Ds);
    if(!fsw)
    Dsw.assign(ns,0.0);
    for(int k=0; k<nz; ++k)
    tc[k] = (double(k)+0.5)/double(nz);

    auto table = [&](double x, double t)
    {
        const double fi = std::max(0.0,std::min((x/L-s[0])/ds,double(ns)-1.000001));
        const int i0 = int(fi);
        const double wi = fi-double(i0);
        for(int k=0; k<nz; ++k)
        {
            row0[k] = dplane[i0*nz+k];
            row1[k] = dplane[(i0+1)*nz+k];
        }
        return (1.0-wi)*interp(t,tc,row0) + wi*interp(t,tc,row1);
    };

    vec xa(nst+ntail+1);
    for(int i=0; i<nst; ++i)
    xa[i] = xmin + (xmax-xmin)*0.5*(1.0-cos(PI*(double(i)+0.5)/double(nst)));
    for(int i=0; i<=ntail; ++i)
    xa[nst+i] = xmax + 2.0*eps + 0.9*L*pow(double(i)/double(ntail),1.6);

    std::vector<curve> cv(xa.size());
    std::vector<std::vector<char> > cw(xa.size());
    double totmax=0.0;
    vec a;
    for(size_t i=0; i<xa.size(); ++i)
    {
        full(xa[i],fb,cv[i],cw[i]);
        arcs_of(cv[i],a);
        totmax = std::max(totmax,a.back());
    }
    vec arcs(nn+1);
    for(int k=0; k<=nn; ++k)
    arcs[k] = double(k)*totmax/double(nn);

    auto station = [&](double x, const curve &pts, const std::vector<char> &wl, curve &q)
    {
        const double total = resample_abs(pts,arcs,q);
        curve wp = pts;
        std::vector<char> ww = wl;
        cut_at_water(wp,ww);
        vec aw;
        arcs_of(wp,aw);
        const double gw = aw.back();
        vec dd(nn+1), ae(nn+1), nry(nn+1,0.0), nrz(nn+1,0.0);
        for(int k=0; k<=nn; ++k)
        {
            ae[k] = std::min(arcs[k],total);
            dd[k] = table(x,std::max(0.0,std::min(ae[k]/std::max(gw,1.0e-12),1.0)));
            const int k0 = std::max(k-1,0), k1 = std::min(k+1,nn);
            const double ty = (q[k1].y-q[k0].y)/double(k1-k0), tz = (q[k1].z-q[k0].z)/double(k1-k0);
            const double ln = sqrt(ty*ty + tz*tz);
            if(ln>1.0e-12)
            {
                nry[k] = tz/ln;                               // outward in the section plane
                nrz[k] = -ty/ln;
            }
        }
        vec sv(s.begin(),s.end());
        const double dsx = interp(x/L,sv,Dsw);
        bool anywall=false;
        for(char c : wl)
        anywall = anywall || c;
        if(anywall)
        {
            const double dtop = table(x,1.0);
            int kw=0;
            for(int k=1; k<=nn; ++k)
            if(fabs(ae[k]-gw)<fabs(ae[kw]-gw))
            kw=k;
            const double w = std::min(0.3*T*nry[kw]*nry[kw] + dsx/(1.25*std::max(dtop,0.002)),0.6*B);
            for(int k=0; k<=nn; ++k)
            if(ae[k]>gw+1.0e-9)
            dd[k] = dtop*std::max(0.0,std::min((1.5-(ae[k]-gw)/std::max(w,1.0e-6))/0.5,1.0));
        }
        else
        for(int k=0; k<=nn; ++k)
        dd[k] += 2.0*dsx/T*std::max(0.0,std::min(1.0-std::max(-q[k].z,0.0)/T,1.0));
        for(int k=0; k<=nn; ++k)
        {
            q[k].y = std::max(q[k].y + dd[k]*nry[k],0.0);
            q[k].z += dd[k]*nrz[k];
        }
    };

    vec xs;
    std::vector<curve> Q;
    curve q;
    for(size_t i=0; i<xa.size(); ++i)
    {
        station(xa[i],cv[i],cw[i],q);
        if(xs.empty())                                            // cap at the bow: the first section on the centre plane
        {
            curve c0 = q;
            for(P2 &v : c0)
            v.y = 0.0;
            xs.push_back(xmin-eps);
            Q.push_back(c0);
        }
        if(int(xs.size())==nst+1)                                 // cap at the stern: the last section, as thin as the wake
        {
            const double ze = z_env(xmax);
            curve c1 = Q.back();
            for(P2 &v : c1)
            {
                const double lim = v.z<=fb ? table(xmax+2.0*eps,std::max(0.0,std::min((v.z-ze)/std::max(-ze,1.0e-9),1.0))) : 0.0;
                v.y = std::min(v.y,lim);
            }
            xs.push_back(xmax+eps);
            Q.push_back(c1);
        }
        xs.push_back(xa[i]);
        Q.push_back(q);
    }

    out.clear();
    double vol=0.0, awet=0.0, vwet=0.0;
    auto add = [&](const double *p, const double *q_, const double *r)
    {
        const double ux=q_[0]-p[0], uy=q_[1]-p[1], uz=q_[2]-p[2], vx=r[0]-p[0], vy=r[1]-p[1], vz=r[2]-p[2];
        const double cx=uy*vz-uz*vy, cy=uz*vx-ux*vz, cz=ux*vy-uy*vx;
        if(sqrt(cx*cx + cy*cy + cz*cz)<=1.0e-13)
        return;
        const double gx=(p[0]+q_[0]+r[0])/3.0, gy=(p[1]+q_[1]+r[1])/3.0, gz=(p[2]+q_[2]+r[2])/3.0;
        vol += (gx*cx + gy*cy + gz*cz)/6.0;
        if(gz<0.0)
        {
            awet += 0.5*sqrt(cx*cx + cy*cy + cz*cz);
            vwet += gy*cy/2.0;
        }
        tri3 t = {p[0]+x_fp,p[1]+y_c,p[2]+z_w,q_[0]+x_fp,q_[1]+y_c,q_[2]+z_w,r[0]+x_fp,r[1]+y_c,r[2]+z_w};
        out.push_back(t);
    };
    for(size_t i=0; i+1<xs.size(); ++i)
    {
        auto P = [&](size_t ii, int k, double sgn, double *v)
        {
            v[0]=xs[ii];
            v[1]=sgn*Q[ii][k].y;
            v[2]=Q[ii][k].z;
        };
        for(int k=0; k<nn; ++k)
        {
            double p00[3],p10[3],p01[3],p11[3],m00[3],m10[3],m01[3],m11[3];
            P(i,k,1.0,p00); P(i+1,k,1.0,p10); P(i,k+1,1.0,p01); P(i+1,k+1,1.0,p11);
            if(std::max(std::max(p00[1],p10[1]),std::max(p01[1],p11[1]))<1.0e-9)
            continue;                                             // on the centre plane: no surface
            P(i,k,-1.0,m00); P(i+1,k,-1.0,m10); P(i,k+1,-1.0,m01); P(i+1,k+1,-1.0,m11);
            add(p00,p11,p10); add(p00,p01,p11);
            add(m00,m10,m11); add(m00,m11,m01);
        }
        double a0[3],a1[3],b0[3],b1[3];
        P(i,nn,1.0,a0); P(i+1,nn,1.0,a1); P(i,nn,-1.0,b0); P(i+1,nn,-1.0,b1);          // top strip (normal up)
        if(std::max(a0[1],a1[1])>1.0e-9)
        {
            add(b0,b1,a1); add(b0,a1,a0);
        }
        P(i,0,1.0,a0); P(i+1,0,1.0,a1); P(i,0,-1.0,b0); P(i+1,0,-1.0,b1);              // bottom strip (normal down)
        if(std::max(a0[1],a1[1])>1.0e-9)
        {
            add(b0,a1,b1); add(b0,a0,a1);
        }
    }
    say("  displacement body: %d triangles, volume %.4f m^3, below the still water level: area %.4f m^2, volume %.4f m^3",
        int(out.size()),vol,awet,vwet);
}

void nhflow_shipbl::run(const std::vector<tri3> &hull, std::vector<tri3> &dispbody)
{
    set_hull(hull);
    build_grid();
    double_body();
    solve();
    body(dispbody);
}
